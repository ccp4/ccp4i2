# Sharing one ACR across multiple instances

By default, each CCP4i2 instance created from `infrastructure.bicep` provisions its own Azure Container Registry. This is simple but wasteful when you run several instances: every instance needs its own image build pipeline, and the base (CCP4 + ARP/wARP) layers are duplicated.

This document describes how to point an instance at a **shared ACR** owned by a different resource group, so that one build pipeline feeds many instances.

> **Current layout (as of 2026-04-20):** `ccp4acrshareduk14fb` in `ccp4i2-shared-rg-uksouth` is the shared ACR. Both `ccp4i2-bicep-*` (DDUDatabase) and `ccp4i2-demo-*` pull from it via identity auth + cross-VNet private endpoints. The legacy `ccp4acrukbwmx` in `ccp4i2-bicep-rg-uksouth` still exists as a safety net but should be deleted after a quiet week.

---

## How cross-instance ACR pulls work

A Container App pulls from `<acrname>.azurecr.io`. If the ACR's public networking is disabled or its private endpoint's DNS zone is linked to the consumer VNet, the pull goes over a private link. For the pull to succeed from a VNet that is not the ACR's home VNet, three things must be in place:

1. **Identity**: the Container App's managed identity must hold `AcrPull` on the ACR resource.
2. **Network path**: the Container App's VNet must have a private endpoint for the shared ACR (cross-VNet peering is also possible but more intrusive).
3. **DNS**: the VNet's `privatelink.azurecr.io` private zone must include A records for the shared ACR pointing at the private endpoint's IP. Without this, the ACR's public CNAME chain (`*.azurecr.io` → `*.privatelink.azurecr.io`) resolves to NXDOMAIN inside the VNet and the pull fails with *"no such host"*.

---

## Wiring a new instance to the shared ACR

Replace the placeholders with your values:

```bash
SHARED_ACR_NAME=ccp4acrshareduk14fb
SHARED_ACR_RG=ccp4i2-shared-rg-uksouth
SHARED_ACR_ID=$(az acr show -g "$SHARED_ACR_RG" -n "$SHARED_ACR_NAME" --query id -o tsv)

CONSUMER_RG=ccp4i2-demo-rg-uksouth
CONSUMER_VNET=ccp4i2-demo-vnet-uk
CONSUMER_PE_SUBNET=private-endpoints-subnet
CONSUMER_MI_PRINCIPAL_ID=bf97f232-06a0-49cf-8c7b-0fc36ddc998c
CONSUMER_MI_ID=/subscriptions/.../ccp4i2-demo-uk-identity
```

### 1. Grant AcrPull to the consumer's managed identity

```bash
az role assignment create \
  --assignee "$CONSUMER_MI_PRINCIPAL_ID" \
  --role AcrPull \
  --scope "$SHARED_ACR_ID"
```

### 2. Create a private endpoint inside the consumer VNet

```bash
az network private-endpoint create \
  --name "${SHARED_ACR_NAME}-pe-from-${CONSUMER_RG}" \
  --resource-group "$CONSUMER_RG" \
  --vnet-name "$CONSUMER_VNET" \
  --subnet "$CONSUMER_PE_SUBNET" \
  --private-connection-resource-id "$SHARED_ACR_ID" \
  --group-id registry \
  --connection-name "${SHARED_ACR_NAME}-conn"
```

### 3. Register the private endpoint in the consumer's `privatelink.azurecr.io` zone

If the consumer already has a `privatelink.azurecr.io` private DNS zone (most do, because they provisioned their own ACR originally):

```bash
az network private-endpoint dns-zone-group create \
  --resource-group "$CONSUMER_RG" \
  --endpoint-name "${SHARED_ACR_NAME}-pe-from-${CONSUMER_RG}" \
  --name default \
  --private-dns-zone privatelink.azurecr.io \
  --zone-name privatelink-azurecr-io
```

This adds two A records to the zone: one for the registry endpoint and one for the regional data endpoint. Verify:

```bash
az network private-dns record-set a list \
  -g "$CONSUMER_RG" \
  --zone-name privatelink.azurecr.io -o table
```

### 4. Point the Container Apps at the shared ACR

For each Container App (web, server, worker):

```bash
# Add the shared ACR as a second registry (identity-based auth, no password)
az containerapp registry set \
  --name <app-name> --resource-group "$CONSUMER_RG" \
  --server "${SHARED_ACR_NAME}.azurecr.io" \
  --identity "$CONSUMER_MI_ID"

# Flip the image reference — this triggers a revision that pulls the new image
az containerapp update \
  --name <app-name> --resource-group "$CONSUMER_RG" \
  --image "${SHARED_ACR_NAME}.azurecr.io/ccp4i2/<web|server>:<tag>"
```

**Important ordering:** `registry set` triggers a revision immediately. If you change the registry first without providing a new image, the new revision tries to pull the *old* image tag from the *new* ACR and fails. Always add the new registry before you remove the old one; flip the image in the same session.

### 5. Drop the old ACR registration

Once the new revisions are Healthy and serving traffic:

```bash
az containerapp registry remove \
  --name <app-name> --resource-group "$CONSUMER_RG" \
  --server <old-acr>.azurecr.io
```

### 6. (Optional) Delete the instance's own ACR

```bash
az acr delete --name <old-acr> --resource-group "$CONSUMER_RG" --yes
```

---

## Known gotchas

- **DNS-only fixes don't work.** Adding a static A record for `<shared-acr>.privatelink.azurecr.io` pointing at the shared ACR's public IP will work today but is brittle: ACR's public IP can change, and the record drifts from reality. Always use a real private endpoint.
- **`registry set` without matching image update fails.** See step 4. The failed revision sticks around on the app until replaced — you'll see "provisioningState: Failed" and the active revision stays on the old image.
- **Admin user / password auth is incompatible with cross-RG ACR.** The default `applications.bicep` uses `containerRegistry.listCredentials()` which only works when the ACR is in the deployment RG. Identity-based auth (this document) is the portable path — we should retire password auth in the bicep.

---

## TODO: bicep integration

The steps in this document are all imperative `az` calls. They should be expressed declaratively in `applications.bicep`:

- Accept `sharedAcrResourceId: string?` parameter (optional).
- If provided, skip the `containerRegistry existing` lookup, skip the password secret, and add the registry stanza with `identity: containerAppsIdentityId`.
- Add a `sharedAcrPrivateEndpoint` resource wiring the PE in the instance's `private-endpoints-subnet`.
- Add a DNS zone group referencing the existing `privatelink.azurecr.io` private zone.
- Emit a role assignment on the consumer MI against the shared ACR.

Once that's done, every new instance created via `deploy-applications.sh` with `--shared-acr-id …` gets the shared-ACR wiring for free, with no imperative post-processing.
