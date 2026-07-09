# macOS code signing + notarisation setup

The release workflow (`.github/workflows/release.yml`) signs and notarises the
macOS app **only when five Apple secrets are present**. If they're absent it
ships an **unsigned** build (testers must run `xattr -cr` — see
[give-it-a-try.md](give-it-a-try.md)). This is a one-time setup that makes every
future `v*` release produce a signed, notarised `.dmg` that opens with no
Gatekeeper warning.

You need an **Apple Developer account** (the CCP4 org's) and **admin** on
`github.com/ccp4/ccp4i2`.

## The five secrets (exact names, case-sensitive)

Add all five as **Repository secrets** on `ccp4/ccp4i2`
(Settings → Secrets and variables → **Actions** → **Repository secrets** →
*New repository secret*).

> **Critical — scope:** they must be **Repository** (or Organization) secrets,
> **not Environment** secrets. The `build-desktop` job runs with no
> `environment:`, so it cannot read secrets scoped to the `pypi`/`npm`
> environments. (This is the mistake that shipped a1–a4 unsigned.)

| Secret | What it is | How to produce the value |
|---|---|---|
| `MAC_CSC_LINK` | Developer ID Application certificate, base64-encoded | see **1** below |
| `MAC_CSC_KEY_PASSWORD` | password you set when exporting the `.p12` | the export password, verbatim |
| `APPLE_API_KEY_P8` | App Store Connect API key, raw `.p8` contents | see **2** below |
| `APPLE_API_KEY_ID` | that key's Key ID (10 chars, e.g. `2X9R4HXF34`) | shown in App Store Connect |
| `APPLE_API_ISSUER` | your team's Issuer ID (a UUID) | shown in App Store Connect |

The gate that flips signed-vs-unsigned checks `MAC_CSC_LINK` **and**
`APPLE_API_KEY_P8`; the other three are consumed during notarisation.

## 1. Signing certificate → `MAC_CSC_LINK` + `MAC_CSC_KEY_PASSWORD`

1. In the Apple Developer portal (or Xcode → Settings → Accounts → Manage
   Certificates), create a **Developer ID Application** certificate. (This is the
   distribution-outside-the-App-Store type — *not* "Apple Development" or "Mac
   App Distribution".)
2. In **Keychain Access** on a Mac that has the cert + its private key: select
   the certificate **and** its private key, right-click → **Export 2 items…** →
   save as a `.p12`. Set an export password — this becomes `MAC_CSC_KEY_PASSWORD`.
3. Base64-encode the `.p12` (electron-builder decodes `CSC_LINK` from base64):
   ```bash
   base64 -i DeveloperIDApplication.p12 | pbcopy   # now paste as MAC_CSC_LINK
   ```
   Paste the base64 blob as `MAC_CSC_LINK`, and the export password as
   `MAC_CSC_KEY_PASSWORD`.

## 2. Notarisation API key → `APPLE_API_KEY_P8` + `_ID` + `_ISSUER`

1. In **App Store Connect** → Users and Access → **Integrations** → **App Store
   Connect API** → **Keys**, create a key with the **Developer** role (sufficient
   for notarisation). Download the `.p8` — **you can only download it once**.
2. From that page:
   - **Key ID** (10 chars) → `APPLE_API_KEY_ID`
   - **Issuer ID** (UUID, top of the Keys page) → `APPLE_API_ISSUER`
3. The **entire contents** of the `.p8` file (including the
   `-----BEGIN PRIVATE KEY-----` / `-----END PRIVATE KEY-----` lines) →
   `APPLE_API_KEY_P8`. The workflow writes it back to a file verbatim, so paste
   it as-is, newlines and all:
   ```bash
   pbcopy < AuthKey_XXXXXXXXXX.p8    # paste as APPLE_API_KEY_P8
   ```

## 3. Verify before cutting a release

The Configure-signing step prints a **presence readout** (names only, never
values). Trigger the release workflow via **workflow_dispatch** (Actions →
Release → Run workflow) and read the macOS job log — you want:

```
Apple signing secrets seen by this job (set = non-empty):
  MAC_CSC_LINK = set (NNNN chars)
  MAC_CSC_KEY_PASSWORD = set (NN chars)
  APPLE_API_KEY_P8 = set (NNN chars)
  APPLE_API_KEY_ID = set (10 chars)
  APPLE_API_ISSUER = set (36 chars)
macOS: signing + notarisation ENABLED
```

Any `EMPTY` line points at a wrong name/scope/repo. Once all five read *set* and
the build logs `signing + notarisation ENABLED` (and later a successful notarise,
not `skipped`), the next `vX.Y.Z[aN]` tag produces a properly signed release.

## Notes

- **Team membership:** the signing cert and the API key must belong to the same
  Apple **team** whose Developer ID will appear on the app.
- **Cert expiry:** Developer ID certs are valid ~5 years; renewing means
  re-exporting a new `.p12` and updating `MAC_CSC_LINK` + `MAC_CSC_KEY_PASSWORD`.
- **Windows:** this covers macOS only. Windows builds are currently unsigned too
  (a separate Authenticode cert would be needed); testers may see a SmartScreen
  prompt. Out of scope here.
