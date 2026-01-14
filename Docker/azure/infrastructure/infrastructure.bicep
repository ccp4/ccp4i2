@description('The location for all resources')
param location string = resourceGroup().location

@description('Resource naming prefix')
param prefix string = 'ccp4i2-bicep'

@description('Environment suffix (e.g., dev, staging, prod)')
param environment string = 'ne'

@description('PostgreSQL administrator password')
@secure()
param postgresAdminPassword string

@description('Skip PostgreSQL deployment if it already exists (use when Azure has transient issues)')
param skipPostgresDeployment bool = false

// Variables
var resourceSuffix = '${prefix}-${environment}'
var acrName = 'ccp4acr${environment}${substring(uniqueString(resourceGroup().id), 0, 4)}' // Keep it short and unique
// CCP4 storage: existing stornekmayz3n2 - PUBLIC access for CCP4 software (working installation with pip packages)
var ccp4StorageAccountName = 'stor${environment}${substring(uniqueString(resourceGroup().id), 0, 8)}'
// Private storage: NEW account for sensitive user data (projects, media, static)
var privateStorageAccountName = 'storprv${environment}${substring(uniqueString(resourceGroup().id), 0, 5)}'
var keyVaultName = 'kv-${environment}-${substring(uniqueString(resourceGroup().id), 0, 6)}' // Keep under 24 chars
var postgresServerName = '${prefix}-db-${environment}'
var containerAppsEnvironmentName = '${prefix}-env-${environment}'
var vnetName = '${prefix}-vnet-${environment}'

// Network Configuration
var vnetAddressSpace = '10.0.0.0/16'
var containerAppsSubnetAddressPrefix = '10.0.0.0/23' // Minimum /23 required for Container Apps
var privateEndpointsSubnetAddressPrefix = '10.0.2.0/24'
var managementSubnetAddressPrefix = '10.0.3.0/24'

// Virtual Network
resource vnet 'Microsoft.Network/virtualNetworks@2023-05-01' = {
  name: vnetName
  location: location
  properties: {
    addressSpace: {
      addressPrefixes: [vnetAddressSpace]
    }
    subnets: [
      {
        name: 'container-apps-subnet'
        properties: {
          addressPrefix: containerAppsSubnetAddressPrefix
          // Delegation will be automatically added by Container Apps environment
        }
      }
      {
        name: 'private-endpoints-subnet'
        properties: {
          addressPrefix: privateEndpointsSubnetAddressPrefix
          privateEndpointNetworkPolicies: 'Disabled'
          privateLinkServiceNetworkPolicies: 'Enabled'
        }
      }
      {
        name: 'management-subnet'
        properties: {
          addressPrefix: managementSubnetAddressPrefix
        }
      }
    ]
  }
}

// Network Security Group for Container Apps
resource containerAppsNsg 'Microsoft.Network/networkSecurityGroups@2023-05-01' = {
  name: '${resourceSuffix}-container-apps-nsg'
  location: location
  properties: {
    securityRules: [
      {
        name: 'AllowHTTPS'
        properties: {
          protocol: 'Tcp'
          sourcePortRange: '*'
          destinationPortRange: '443'
          sourceAddressPrefix: '*'
          destinationAddressPrefix: '*'
          access: 'Allow'
          priority: 1000
          direction: 'Inbound'
        }
      }
      {
        name: 'AllowHTTP'
        properties: {
          protocol: 'Tcp'
          sourcePortRange: '*'
          destinationPortRange: '80'
          sourceAddressPrefix: '*'
          destinationAddressPrefix: '*'
          access: 'Allow'
          priority: 1010
          direction: 'Inbound'
        }
      }
    ]
  }
}

// Network Security Group for Private Endpoints
resource privateEndpointsNsg 'Microsoft.Network/networkSecurityGroups@2023-05-01' = {
  name: '${resourceSuffix}-private-endpoints-nsg'
  location: location
  properties: {
    securityRules: [
      {
        name: 'DenyAll'
        properties: {
          protocol: '*'
          sourcePortRange: '*'
          destinationPortRange: '*'
          sourceAddressPrefix: '*'
          destinationAddressPrefix: '*'
          access: 'Deny'
          priority: 4096
          direction: 'Inbound'
        }
      }
    ]
  }
}

// Private DNS Zones
resource privateDnsZoneStorage 'Microsoft.Network/privateDnsZones@2020-06-01' = {
  name: 'privatelink.file.${az.environment().suffixes.storage}'
  location: 'global'
}

resource privateDnsZoneBlob 'Microsoft.Network/privateDnsZones@2020-06-01' = {
  name: 'privatelink.blob.${az.environment().suffixes.storage}'
  location: 'global'
}

resource privateDnsZoneKeyVault 'Microsoft.Network/privateDnsZones@2020-06-01' = {
  name: 'privatelink.vaultcore.azure.net'
  location: 'global'
}

resource privateDnsZoneAcr 'Microsoft.Network/privateDnsZones@2020-06-01' = {
  name: 'privatelink.azurecr.io'
  location: 'global'
}

resource privateDnsZonePostgres 'Microsoft.Network/privateDnsZones@2020-06-01' = {
  name: 'privatelink.postgres.database.azure.com'
  location: 'global'
}

resource privateDnsZoneServiceBus 'Microsoft.Network/privateDnsZones@2020-06-01' = {
  name: 'privatelink.servicebus.windows.net'
  location: 'global'
}

// Link Private DNS Zones to VNet
resource privateDnsZoneStorageLink 'Microsoft.Network/privateDnsZones/virtualNetworkLinks@2020-06-01' = {
  name: 'storage-link'
  parent: privateDnsZoneStorage
  location: 'global'
  properties: {
    registrationEnabled: false
    virtualNetwork: {
      id: vnet.id
    }
  }
}

resource privateDnsZoneBlobLink 'Microsoft.Network/privateDnsZones/virtualNetworkLinks@2020-06-01' = {
  name: 'blob-link'
  parent: privateDnsZoneBlob
  location: 'global'
  properties: {
    registrationEnabled: false
    virtualNetwork: {
      id: vnet.id
    }
  }
}

resource privateDnsZoneKeyVaultLink 'Microsoft.Network/privateDnsZones/virtualNetworkLinks@2020-06-01' = {
  name: 'keyvault-link'
  parent: privateDnsZoneKeyVault
  location: 'global'
  properties: {
    registrationEnabled: false
    virtualNetwork: {
      id: vnet.id
    }
  }
}

resource privateDnsZoneAcrLink 'Microsoft.Network/privateDnsZones/virtualNetworkLinks@2020-06-01' = {
  name: 'acr-link'
  parent: privateDnsZoneAcr
  location: 'global'
  properties: {
    registrationEnabled: false
    virtualNetwork: {
      id: vnet.id
    }
  }
}

resource privateDnsZonePostgresLink 'Microsoft.Network/privateDnsZones/virtualNetworkLinks@2020-06-01' = {
  name: 'postgres-link'
  parent: privateDnsZonePostgres
  location: 'global'
  properties: {
    registrationEnabled: false
    virtualNetwork: {
      id: vnet.id
    }
  }
}

resource privateDnsZoneServiceBusLink 'Microsoft.Network/privateDnsZones/virtualNetworkLinks@2020-06-01' = {
  name: 'servicebus-link'
  parent: privateDnsZoneServiceBus
  location: 'global'
  properties: {
    registrationEnabled: false
    virtualNetwork: {
      id: vnet.id
    }
  }
}

// Existing resources (deployed separately)
resource serviceBusNamespace 'Microsoft.ServiceBus/namespaces@2024-01-01' existing = {
  name: '${prefix}-servicebus'
}
resource containerRegistry 'Microsoft.ContainerRegistry/registries@2023-07-01' = {
  name: acrName
  location: location
  sku: {
    name: 'Premium' // Premium SKU required for private endpoints
  }
  properties: {
    adminUserEnabled: true
    publicNetworkAccess: 'Enabled' // Start with public access, can be disabled later
    networkRuleBypassOptions: 'AzureServices' // Allow Azure services like ACR Tasks
    networkRuleSet: {
      defaultAction: 'Allow' // Allow all access when public is enabled
    }
  }
}

// Storage Account for CCP4 software (stornekmayz3n2 - existing working installation)
// PUBLIC access - contains read-only CCP4 installation with pip packages, not sensitive user data
resource ccp4StorageAccount 'Microsoft.Storage/storageAccounts@2023-01-01' = {
  name: ccp4StorageAccountName
  location: location
  sku: {
    name: 'Standard_LRS'
  }
  kind: 'StorageV2'
  properties: {
    allowBlobPublicAccess: false  // No anonymous blob access, SAS/key required
    allowSharedKeyAccess: true
    minimumTlsVersion: 'TLS1_2'
    supportsHttpsTrafficOnly: true
    publicNetworkAccess: 'Enabled'
    networkAcls: {
      defaultAction: 'Allow'  // Allow access from anywhere with valid SAS/key
      bypass: 'AzureServices'
    }
  }
}

// Storage Account (PRIVATE) - for sensitive user data (projects, media, static)
// Network access is restricted - only accessible via:
// - Private endpoints (container apps within VNet)
// - Azure services bypass
// - Temporarily added IPs (via admin scripts)
resource privateStorageAccount 'Microsoft.Storage/storageAccounts@2023-01-01' = {
  name: privateStorageAccountName
  location: location
  sku: {
    name: 'Standard_LRS'
  }
  kind: 'StorageV2'
  properties: {
    allowBlobPublicAccess: false  // No anonymous access
    allowSharedKeyAccess: true
    minimumTlsVersion: 'TLS1_2'
    supportsHttpsTrafficOnly: true
    publicNetworkAccess: 'Enabled'  // Required for IP-based firewall rules to work
    networkAcls: {
      defaultAction: 'Deny'  // Deny by default - only private endpoints and allowed IPs
      bypass: 'AzureServices'
      ipRules: []  // Admin scripts will temporarily add IPs as needed
      virtualNetworkRules: []
    }
  }
}

// File Share for CCP4 software (on CCP4 storage account - existing stornekmayz3n2)
resource ccp4dataShare 'Microsoft.Storage/storageAccounts/fileServices/shares@2023-01-01' = {
  name: '${ccp4StorageAccount.name}/default/ccp4data'
  properties: {
    shareQuota: 100
  }
}

// File shares on PRIVATE storage account (for sensitive user data)
resource staticfilesShare 'Microsoft.Storage/storageAccounts/fileServices/shares@2023-01-01' = {
  name: '${privateStorageAccount.name}/default/staticfiles'
  properties: {
    shareQuota: 10
  }
}

resource mediafilesShare 'Microsoft.Storage/storageAccounts/fileServices/shares@2023-01-01' = {
  name: '${privateStorageAccount.name}/default/mediafiles'
  properties: {
    shareQuota: 50
  }
}

// File share for CCP4i2 projects (on private storage - separate from ccp4data)
// Uses Azure Files for POSIX filesystem semantics needed by CCP4 tools
resource projectsShare 'Microsoft.Storage/storageAccounts/fileServices/shares@2023-01-01' = {
  name: '${privateStorageAccount.name}/default/ccp4i2-projects'
  properties: {
    shareQuota: 5120  // 5 TiB quota for projects (matching legacy storage)
  }
}

// Blob services on PRIVATE storage account
resource blobServices 'Microsoft.Storage/storageAccounts/blobServices@2023-01-01' = {
  name: 'default'
  parent: privateStorageAccount
}

resource stagingUploadsContainer 'Microsoft.Storage/storageAccounts/blobServices/containers@2023-01-01' = {
  name: 'staging-uploads'
  parent: blobServices
  properties: {
    publicAccess: 'None'  // No anonymous access, SAS required
  }
}

// Blob container for Django file uploads (using django-storages)
// Uses Blob Storage (cheaper than Files) for object storage of Django uploads
resource djangoUploadsContainer 'Microsoft.Storage/storageAccounts/blobServices/containers@2023-01-01' = {
  name: 'django-uploads'
  parent: blobServices
  properties: {
    publicAccess: 'None'  // No anonymous access, managed identity auth
  }
}

// Key Vault
resource keyVault 'Microsoft.KeyVault/vaults@2023-02-01' = {
  name: keyVaultName
  location: location
  properties: {
    tenantId: subscription().tenantId
    sku: {
      family: 'A'
      name: 'standard'
    }
    enabledForTemplateDeployment: true
    enableRbacAuthorization: true
    publicNetworkAccess: 'Disabled'
    networkAcls: {
      defaultAction: 'Deny'
      bypass: 'AzureServices'
    }
  }
}

// Store PostgreSQL password in Key Vault
resource postgresPasswordSecret 'Microsoft.KeyVault/vaults/secrets@2023-02-01' = {
  name: 'database-admin-password'
  parent: keyVault
  properties: {
    value: postgresAdminPassword
  }
}

// PostgreSQL Flexible Server - conditionally deployed
// When skipPostgresDeployment is true, reference existing server instead of deploying
resource postgresServer 'Microsoft.DBforPostgreSQL/flexibleServers@2023-03-01-preview' = if (!skipPostgresDeployment) {
  name: postgresServerName
  location: location
  sku: {
    name: 'Standard_B1ms'
    tier: 'Burstable'
  }
  properties: {
    administratorLogin: 'ccp4i2'
    administratorLoginPassword: postgresAdminPassword
    version: '15'
    storage: {
      storageSizeGB: 32
    }
    backup: {
      backupRetentionDays: 7
      geoRedundantBackup: 'Disabled'
    }
    // Network configuration will be handled by private endpoint
  }
}

// Reference existing PostgreSQL server when skipping deployment
resource existingPostgresServer 'Microsoft.DBforPostgreSQL/flexibleServers@2023-03-01-preview' existing = if (skipPostgresDeployment) {
  name: postgresServerName
}

// Remove the old firewall rule since we're using private endpoints

// Private Endpoints for NEW private storage account
// Using new names because existing endpoints can't change their target resource
resource storagePrivateEndpoint 'Microsoft.Network/privateEndpoints@2023-05-01' = {
  name: '${resourceSuffix}-storage-private-pe'
  location: location
  properties: {
    subnet: {
      id: '${vnet.id}/subnets/private-endpoints-subnet'
    }
    privateLinkServiceConnections: [
      {
        name: 'storage-connection'
        properties: {
          privateLinkServiceId: privateStorageAccount.id
          groupIds: ['file']
        }
      }
    ]
  }
}

resource storagePrivateDnsZoneGroup 'Microsoft.Network/privateEndpoints/privateDnsZoneGroups@2023-05-01' = {
  name: 'storage-dns-zone-group'
  parent: storagePrivateEndpoint
  properties: {
    privateDnsZoneConfigs: [
      {
        name: 'storage-config'
        properties: {
          privateDnsZoneId: privateDnsZoneStorage.id
        }
      }
    ]
  }
}

// Private endpoint for Blob storage (django-storages uploads)
resource blobPrivateEndpoint 'Microsoft.Network/privateEndpoints@2023-05-01' = {
  name: '${resourceSuffix}-blob-private-pe'
  location: location
  properties: {
    subnet: {
      id: '${vnet.id}/subnets/private-endpoints-subnet'
    }
    privateLinkServiceConnections: [
      {
        name: 'blob-connection'
        properties: {
          privateLinkServiceId: privateStorageAccount.id
          groupIds: ['blob']
        }
      }
    ]
  }
}

resource blobPrivateDnsZoneGroup 'Microsoft.Network/privateEndpoints/privateDnsZoneGroups@2023-05-01' = {
  name: 'blob-dns-zone-group'
  parent: blobPrivateEndpoint
  properties: {
    privateDnsZoneConfigs: [
      {
        name: 'blob-config'
        properties: {
          privateDnsZoneId: privateDnsZoneBlob.id
        }
      }
    ]
  }
}

resource keyVaultPrivateEndpoint 'Microsoft.Network/privateEndpoints@2023-05-01' = {
  name: '${resourceSuffix}-keyvault-pe'
  location: location
  properties: {
    subnet: {
      id: '${vnet.id}/subnets/private-endpoints-subnet'
    }
    privateLinkServiceConnections: [
      {
        name: 'keyvault-connection'
        properties: {
          privateLinkServiceId: keyVault.id
          groupIds: ['vault']
        }
      }
    ]
  }
}

resource keyVaultPrivateDnsZoneGroup 'Microsoft.Network/privateEndpoints/privateDnsZoneGroups@2023-05-01' = {
  name: 'keyvault-dns-zone-group'
  parent: keyVaultPrivateEndpoint
  properties: {
    privateDnsZoneConfigs: [
      {
        name: 'keyvault-config'
        properties: {
          privateDnsZoneId: privateDnsZoneKeyVault.id
        }
      }
    ]
  }
}

resource acrPrivateEndpoint 'Microsoft.Network/privateEndpoints@2023-05-01' = {
  name: '${resourceSuffix}-acr-pe'
  location: location
  properties: {
    subnet: {
      id: '${vnet.id}/subnets/private-endpoints-subnet'
    }
    privateLinkServiceConnections: [
      {
        name: 'acr-connection'
        properties: {
          privateLinkServiceId: containerRegistry.id
          groupIds: ['registry']
        }
      }
    ]
  }
}

resource acrPrivateDnsZoneGroup 'Microsoft.Network/privateEndpoints/privateDnsZoneGroups@2023-05-01' = {
  name: 'acr-dns-zone-group'
  parent: acrPrivateEndpoint
  properties: {
    privateDnsZoneConfigs: [
      {
        name: 'acr-config'
        properties: {
          privateDnsZoneId: privateDnsZoneAcr.id
        }
      }
    ]
  }
}

resource postgresPrivateEndpoint 'Microsoft.Network/privateEndpoints@2023-05-01' = {
  name: '${resourceSuffix}-postgres-pe'
  location: location
  properties: {
    subnet: {
      id: '${vnet.id}/subnets/private-endpoints-subnet'
    }
    privateLinkServiceConnections: [
      {
        name: 'postgres-connection'
        properties: {
          privateLinkServiceId: skipPostgresDeployment ? existingPostgresServer.id : postgresServer.id
          groupIds: ['postgresqlServer']
        }
      }
    ]
  }
}

resource postgresPrivateDnsZoneGroup 'Microsoft.Network/privateEndpoints/privateDnsZoneGroups@2023-05-01' = {
  name: 'postgres-dns-zone-group'
  parent: postgresPrivateEndpoint
  properties: {
    privateDnsZoneConfigs: [
      {
        name: 'postgres-config'
        properties: {
          privateDnsZoneId: privateDnsZonePostgres.id
        }
      }
    ]
  }
}

resource serviceBusPrivateEndpoint 'Microsoft.Network/privateEndpoints@2023-05-01' = {
  name: '${resourceSuffix}-servicebus-pe'
  location: location
  properties: {
    subnet: {
      id: '${vnet.id}/subnets/private-endpoints-subnet'
    }
    privateLinkServiceConnections: [
      {
        name: 'servicebus-connection'
        properties: {
          privateLinkServiceId: serviceBusNamespace.id
          groupIds: ['namespace']
        }
      }
    ]
  }
}

resource serviceBusPrivateDnsZoneGroup 'Microsoft.Network/privateEndpoints/privateDnsZoneGroups@2023-05-01' = {
  name: 'servicebus-dns-zone-group'
  parent: serviceBusPrivateEndpoint
  properties: {
    privateDnsZoneConfigs: [
      {
        name: 'servicebus-config'
        properties: {
          privateDnsZoneId: privateDnsZoneServiceBus.id
        }
      }
    ]
  }
}

// Container Apps Environment
resource containerAppsEnvironment 'Microsoft.App/managedEnvironments@2023-05-01' = {
  name: containerAppsEnvironmentName
  location: location
  properties: {
    vnetConfiguration: {
      infrastructureSubnetId: '${vnet.id}/subnets/container-apps-subnet'
      internal: false // Set to true for completely internal environment
    }
    appLogsConfiguration: {
      destination: 'log-analytics'
      logAnalyticsConfiguration: {
        customerId: logAnalyticsWorkspace.properties.customerId
        sharedKey: logAnalyticsWorkspace.listKeys().primarySharedKey
      }
    }
  }
}

// Log Analytics Workspace for Container Apps
resource logAnalyticsWorkspace 'Microsoft.OperationalInsights/workspaces@2022-10-01' = {
  name: '${resourceSuffix}-logs'
  location: location
  properties: {
    sku: {
      name: 'PerGB2018'
    }
    retentionInDays: 30
  }
}

// Shared User-Assigned Managed Identity for Container Apps
resource containerAppsIdentity 'Microsoft.ManagedIdentity/userAssignedIdentities@2023-01-31' = {
  name: '${resourceSuffix}-identity'
  location: location
}

// Role Assignments for Container Apps Identity on Private Storage Account
// Storage Blob Data Contributor - allows read/write access to blobs
resource storageBlobDataContributorRole 'Microsoft.Authorization/roleAssignments@2022-04-01' = {
  name: guid(privateStorageAccount.id, containerAppsIdentity.id, 'ba92f5b4-2d11-453d-a403-e96b0029c9fe')
  scope: privateStorageAccount
  properties: {
    roleDefinitionId: subscriptionResourceId('Microsoft.Authorization/roleDefinitions', 'ba92f5b4-2d11-453d-a403-e96b0029c9fe')
    principalId: containerAppsIdentity.properties.principalId
    principalType: 'ServicePrincipal'
  }
}

// Storage Blob Delegator - allows generating user delegation SAS tokens
resource storageBlobDelegatorRole 'Microsoft.Authorization/roleAssignments@2022-04-01' = {
  name: guid(privateStorageAccount.id, containerAppsIdentity.id, 'db58b8e5-c6ad-4a2a-8342-4190687cbf4a')
  scope: privateStorageAccount
  properties: {
    roleDefinitionId: subscriptionResourceId('Microsoft.Authorization/roleDefinitions', 'db58b8e5-c6ad-4a2a-8342-4190687cbf4a')
    principalId: containerAppsIdentity.properties.principalId
    principalType: 'ServicePrincipal'
  }
}

// Storage for Container Apps Environment - CCP4 software (from existing stornekmayz3n2)
// Using new mount name because existing mount can't change storage account properties
resource containerAppsCcp4Storage 'Microsoft.App/managedEnvironments/storages@2023-05-01' = {
  name: 'ccp4-software'
  parent: containerAppsEnvironment
  properties: {
    azureFile: {
      accountName: ccp4StorageAccount.name
      accountKey: ccp4StorageAccount.listKeys().keys[0].value
      shareName: 'ccp4data'
      accessMode: 'ReadOnly'  // CCP4 software is read-only
    }
  }
}

// Storage mounts for private user data (from new secure storage account)
// Using new mount names with '-private' suffix because existing mounts can't change storage accounts
resource containerAppsStaticStorage 'Microsoft.App/managedEnvironments/storages@2023-05-01' = {
  name: 'staticfiles-private'
  parent: containerAppsEnvironment
  properties: {
    azureFile: {
      accountName: privateStorageAccount.name
      accountKey: privateStorageAccount.listKeys().keys[0].value
      shareName: 'staticfiles'
      accessMode: 'ReadWrite'
    }
  }
}

resource containerAppsMediaStorage 'Microsoft.App/managedEnvironments/storages@2023-05-01' = {
  name: 'mediafiles-private'
  parent: containerAppsEnvironment
  properties: {
    azureFile: {
      accountName: privateStorageAccount.name
      accountKey: privateStorageAccount.listKeys().keys[0].value
      shareName: 'mediafiles'
      accessMode: 'ReadWrite'
    }
  }
}

// Storage mount for CCP4i2 projects (on private storage - separate from ccp4data)
resource containerAppsProjectsStorage 'Microsoft.App/managedEnvironments/storages@2023-05-01' = {
  name: 'projects-private'
  parent: containerAppsEnvironment
  properties: {
    azureFile: {
      accountName: privateStorageAccount.name
      accountKey: privateStorageAccount.listKeys().keys[0].value
      shareName: 'ccp4i2-projects'
      accessMode: 'ReadWrite'
    }
  }
}

// Outputs
output acrName string = containerRegistry.name
output acrLoginServer string = containerRegistry.properties.loginServer
output ccp4StorageAccountName string = ccp4StorageAccount.name
output privateStorageAccountName string = privateStorageAccount.name
output keyVaultName string = keyVault.name
output postgresServerName string = skipPostgresDeployment ? existingPostgresServer.name : postgresServer.name
output postgresServerFqdn string = skipPostgresDeployment ? existingPostgresServer.properties.fullyQualifiedDomainName : postgresServer.properties.fullyQualifiedDomainName
output containerAppsEnvironmentName string = containerAppsEnvironment.name
output containerAppsEnvironmentId string = containerAppsEnvironment.id
output vnetId string = vnet.id
output vnetName string = vnet.name
output containerAppsSubnetId string = '${vnet.id}/subnets/container-apps-subnet'
output privateEndpointsSubnetId string = '${vnet.id}/subnets/private-endpoints-subnet'
output resourceGroupName string = resourceGroup().name
output containerAppsIdentityId string = containerAppsIdentity.id
output containerAppsIdentityPrincipalId string = containerAppsIdentity.properties.principalId
output containerAppsIdentityClientId string = containerAppsIdentity.properties.clientId
