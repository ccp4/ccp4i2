@description('Container Apps Environment ID')
param containerAppsEnvironmentId string

@description('Azure Container Registry login server')
param acrLoginServer string

@description('Azure Container Registry name')
param acrName string

@description('PostgreSQL server FQDN - will resolve to private IP via private DNS zone')
param postgresServerFqdn string

@description('Key Vault name - accessed via private endpoint')
param keyVaultName string

@description('Web container image tag')
param imageTagWeb string = 'latest'

@description('Server container image tag')
param imageTagServer string = 'latest'

@description('Resource naming prefix')
param prefix string = 'ccp4i2-bicep'

@description('Azure AD Client ID for frontend authentication')
param aadClientId string = '386da83f-1bf4-4ad8-b742-79b600e2208b'

@description('Azure AD Tenant ID for frontend authentication')
param aadTenantId string = '9c5012c9-b616-44c2-a917-66814fbe3e87'

@description('Bootstrap platform admin emails (comma-separated)')
param platformAdminEmails string = ''

@description('Shared Container Apps Identity ID')
param containerAppsIdentityId string

@description('Shared Container Apps Identity Principal ID (for RBAC)')
param containerAppsIdentityPrincipalId string

@description('Shared Container Apps Identity Client ID (for DefaultAzureCredential)')
param containerAppsIdentityClientId string

@description('CCP4 version directory name (e.g., ccp4-9, ccp4-20251105)')
param ccp4Version string = 'ccp4-20251105'

@description('Storage account name for staged uploads (SAS URL generation)')
param storageAccountName string

@description('Skip CCP4 storage mount if CCP4 is baked into container image')
param skipCcp4Storage bool = false

// - PostgreSQL is accessed via private endpoint (no public access)
// - Key Vault is accessed via private endpoint (no public access)
// - Storage Account is accessed via private endpoint (no public access)
// - Container Apps run in dedicated VNet subnet with proper delegation
// - All communication happens within the Azure private network

// Variables
var serverAppName = '${prefix}-server'
var webAppName = '${prefix}-web'
var workerAppName = '${prefix}-worker'

// Conditional CCP4 volume configuration
// When skipCcp4Storage=true, CCP4 is baked into the container image at /mnt/ccp4data
var ccp4Volume = skipCcp4Storage ? [] : [
  {
    name: 'ccp4data-volume'
    storageName: 'ccp4-software'
    storageType: 'AzureFile'
  }
]
var ccp4VolumeMount = skipCcp4Storage ? [] : [
  {
    volumeName: 'ccp4data-volume'
    mountPath: '/mnt/ccp4data'
  }
]

// Existing resources
resource containerRegistry 'Microsoft.ContainerRegistry/registries@2023-07-01' existing = {
  name: acrName
}

resource keyVault 'Microsoft.KeyVault/vaults@2023-02-01' existing = {
  name: keyVaultName
}

// Reference existing Container Apps Environment to get its default domain
resource containerAppsEnvironment 'Microsoft.App/managedEnvironments@2023-05-01' existing = {
  name: last(split(containerAppsEnvironmentId, '/'))
}

// Internal URL for server (accessible only within the Container Apps environment)
var serverInternalFqdn = '${serverAppName}.internal.${containerAppsEnvironment.properties.defaultDomain}'

// Dynamic CORS origin using the Container Apps Environment's default domain
// This avoids hardcoding region-specific URLs like 'whitecliff-258bc831.northeurope.azurecontainerapps.io'
var webAppExternalFqdn = '${webAppName}.${containerAppsEnvironment.properties.defaultDomain}'

// Server Container App
resource serverApp 'Microsoft.App/containerApps@2023-05-01' = {
  name: serverAppName
  location: resourceGroup().location
  identity: {
    type: 'UserAssigned'
    userAssignedIdentities: {
      '${containerAppsIdentityId}': {}
    }
  }
  properties: {
    managedEnvironmentId: containerAppsEnvironmentId
    configuration: {
      activeRevisionsMode: 'Single'
      ingress: {
        external: false // Internal only - accessible via web proxy, not directly from internet
        targetPort: 8000
        allowInsecure: true // Allow HTTP for internal VNet communication
        traffic: [
          {
            weight: 100
            latestRevision: true
          }
        ]
      }
      registries: [
        {
          server: acrLoginServer
          username: acrName
          passwordSecretRef: 'registry-password'
        }
      ]
      secrets: [
        {
          name: 'registry-password'
          value: containerRegistry.listCredentials().passwords[0].value
        }
        {
          name: 'db-password'
          keyVaultUrl: '${keyVault.properties.vaultUri}secrets/database-admin-password'
          identity: containerAppsIdentityId
        }
        {
          name: 'django-secret-key'
          keyVaultUrl: '${keyVault.properties.vaultUri}secrets/django-secret-key'
          identity: containerAppsIdentityId
        }
        {
          name: 'servicebus-connection'
          keyVaultUrl: '${keyVault.properties.vaultUri}secrets/servicebus-connection'
          identity: containerAppsIdentityId
        }
      ]
    }
    template: {
      containers: [
        {
          name: 'server'
          image: '${acrLoginServer}/ccp4i2/server:${imageTagServer}'
          command: ['/bin/bash']
          args: ['-c', 'export PYTHONPATH="/mnt/ccp4data/py-packages-${ccp4Version}:/usr/src/app:$PYTHONPATH" && exec /usr/src/app/startup.sh']
          resources: {
            cpu: json('2.0')
            memory: '4.0Gi'
          }
          probes: [
            {
              type: 'Liveness'
              httpGet: {
                path: '/api/ccp4i2/health/'
                port: 8000
              }
              initialDelaySeconds: 60 // Reduced from 120 to comply with Azure limits (max 60)
              periodSeconds: 60 // Check less frequently
              failureThreshold: 3
              successThreshold: 1
              timeoutSeconds: 10 // Add timeout for slow responses
            }
            {
              type: 'Readiness'
              httpGet: {
                path: '/api/ccp4i2/health/'  // Changed from /projects/ which requires auth
                port: 8000
              }
              initialDelaySeconds: 60
              periodSeconds: 30
              failureThreshold: 3
              successThreshold: 1
              timeoutSeconds: 10
            }
            {
              type: 'Startup'
              httpGet: {
                path: '/api/ccp4i2/health/'
                port: 8000
              }
              initialDelaySeconds: 60 // Reduced from 180 to comply with Azure limits (max 60) - Django should start within 60 seconds
              periodSeconds: 30
              failureThreshold: 10 // Allow more failures during startup
              successThreshold: 1
              timeoutSeconds: 10
            }
          ]
          env: [
            {
              name: 'EXECUTION_MODE'
              value: 'azure'
            }
            {
              name: 'DJANGO_SETTINGS_MODULE'
              value: 'azure_extensions.settings'
            }
            {
              name: 'CCP4_VERSION'
              value: ccp4Version
            }
            {
              name: 'DB_HOST'
              value: postgresServerFqdn // Will resolve to private IP via private DNS zone
            }
            {
              name: 'DB_PORT'
              value: '5432'
            }
            {
              name: 'DB_USER'
              value: 'ccp4i2'
            }
            {
              name: 'DB_NAME'
              value: 'postgres'
            }
            {
              name: 'DB_PASSWORD'
              secretRef: 'db-password'
            }
            {
              name: 'SECRET_KEY'
              secretRef: 'django-secret-key'
            }
            {
              name: 'DB_SSL_MODE'
              value: 'require'
            }
            {
              name: 'DB_SSL_ROOT_CERT'
              value: 'true'
            }
            {
              name: 'DB_SSL_REQUIRE_CERT'
              value: 'false' // Private endpoint uses Azure's trusted certificates
            }
            {
              name: 'DEBUG'
              value: 'false'
            }
            {
              name: 'CCP4_DATA_PATH'
              value: '/mnt/ccp4data'
            }
            {
              name: 'CCP4I2_PROJECTS_DIR'
              value: '/mnt/projects'
            }
            {
              name: 'ALLOWED_HOSTS'
              value: '*'  // Allow all hosts - health probes use various internal IPs/hostnames
            }
            {
              name: 'CORS_ALLOWED_ORIGINS'
              value: 'http://${webAppName},https://${webAppExternalFqdn}'
            }
            {
              name: 'CORS_ALLOW_CREDENTIALS'
              value: 'True'
            }
            {
              name: 'SERVICE_BUS_CONNECTION_STRING'
              secretRef: 'servicebus-connection'
            }
            {
              name: 'SERVICE_BUS_QUEUE_NAME'
              value: '${prefix}-jobs'
            }
            {
              name: 'FILE_UPLOAD_MAX_MEMORY_SIZE'
              value: '104857600' // 100MB in bytes
            }
            {
              name: 'DATA_UPLOAD_MAX_MEMORY_SIZE'
              value: '104857600' // 100MB in bytes
            }
            {
              name: 'FILE_UPLOAD_MAX_NUMBER_FILES'
              value: '10'
            }
            // Azure AD authentication configuration
            // When CCP4I2_REQUIRE_AUTH=true, Django will validate JWT tokens
            {
              name: 'CCP4I2_REQUIRE_AUTH'
              value: 'true'
            }
            {
              name: 'AZURE_AD_TENANT_ID'
              value: aadTenantId
            }
            {
              name: 'AZURE_AD_CLIENT_ID'
              value: aadClientId
            }
            {
              name: 'PLATFORM_ADMIN_EMAILS'
              value: platformAdminEmails
            }
            // Azure Storage for staged uploads (large file uploads via SAS URL)
            {
              name: 'AZURE_STORAGE_ACCOUNT_NAME'
              value: storageAccountName
            }
            // Managed Identity Client ID for DefaultAzureCredential
            // Required for User-Assigned Managed Identity to authenticate to Azure services
            {
              name: 'AZURE_CLIENT_ID'
              value: containerAppsIdentityClientId
            }
          ]
          volumeMounts: concat(ccp4VolumeMount, [
            {
              volumeName: 'staticfiles-volume'
              mountPath: '/mnt/staticfiles'
            }
            {
              volumeName: 'mediafiles-volume'
              mountPath: '/mnt/mediafiles'
            }
            {
              volumeName: 'projects-volume'
              mountPath: '/mnt/projects'
            }
          ])
        }
      ]
      scale: {
        minReplicas: 1
        maxReplicas: 10
        rules: [
          {
            name: 'cpu-scaling'
            custom: {
              type: 'cpu'
              metadata: {
                type: 'Utilization'
                value: '70'
              }
            }
          }
          {
            name: 'memory-scaling'
            custom: {
              type: 'memory'
              metadata: {
                type: 'Utilization'
                value: '70'
              }
            }
          }
          {
            name: 'http-scaling'
            http: {
              metadata: {
                concurrentRequests: '20' // Lower threshold for faster scaling on file uploads
              }
            }
          }
        ]
      }
      volumes: concat(ccp4Volume, [
        {
          name: 'staticfiles-volume'
          storageName: 'staticfiles-private'
          storageType: 'AzureFile'
        }
        {
          name: 'mediafiles-volume'
          storageName: 'mediafiles-private'
          storageType: 'AzureFile'
        }
        {
          name: 'projects-volume'
          storageName: 'projects-private'
          storageType: 'AzureFile'
        }
      ])
    }
  }
}

// Worker Container App
resource workerApp 'Microsoft.App/containerApps@2023-05-01' = {
  name: workerAppName
  location: resourceGroup().location
  identity: {
    type: 'UserAssigned'
    userAssignedIdentities: {
      '${containerAppsIdentityId}': {}
    }
  }
  properties: {
    managedEnvironmentId: containerAppsEnvironmentId
    configuration: {
      activeRevisionsMode: 'Single'
      // No ingress - worker is a background job processor, doesn't serve HTTP
      registries: [
        {
          server: acrLoginServer
          username: acrName
          passwordSecretRef: 'registry-password'
        }
      ]
      secrets: [
        {
          name: 'registry-password'
          value: containerRegistry.listCredentials().passwords[0].value
        }
        {
          name: 'db-password'
          keyVaultUrl: '${keyVault.properties.vaultUri}secrets/database-admin-password'
          identity: containerAppsIdentityId
        }
        {
          name: 'django-secret-key'
          keyVaultUrl: '${keyVault.properties.vaultUri}secrets/django-secret-key'
          identity: containerAppsIdentityId
        }
        {
          name: 'servicebus-connection'
          keyVaultUrl: '${keyVault.properties.vaultUri}secrets/servicebus-connection'
          identity: containerAppsIdentityId
        }
      ]
    }
    template: {
      containers: [
        {
          name: 'worker'
          image: '${acrLoginServer}/ccp4i2/server:${imageTagServer}'
          command: ['/bin/bash']
          args: ['-c', 'export PYTHONPATH="/mnt/ccp4data/py-packages-${ccp4Version}:/usr/src/app:$PYTHONPATH" && exec /usr/src/app/startup-worker.sh']
          resources: {
            cpu: json('2.0')  // Maximum for Consumption plan
            memory: '4.0Gi'   // Maximum for Consumption plan (2 vCPU supports up to 4GB)
          }
          env: [
            {
              name: 'EXECUTION_MODE'
              value: 'azure'
            }
            {
              name: 'DJANGO_SETTINGS_MODULE'
              value: 'azure_extensions.settings'
            }
            {
              name: 'CCP4_VERSION'
              value: ccp4Version
            }
            {
              name: 'DB_HOST'
              value: postgresServerFqdn
            }
            {
              name: 'DB_PORT'
              value: '5432'
            }
            {
              name: 'DB_USER'
              value: 'ccp4i2'
            }
            {
              name: 'DB_NAME'
              value: 'postgres'
            }
            {
              name: 'DB_PASSWORD'
              secretRef: 'db-password'
            }
            {
              name: 'SECRET_KEY'
              secretRef: 'django-secret-key'
            }
            {
              name: 'DB_SSL_MODE'
              value: 'require'
            }
            {
              name: 'DB_SSL_ROOT_CERT'
              value: 'true'
            }
            {
              name: 'DB_SSL_REQUIRE_CERT'
              value: 'false'
            }
            {
              name: 'CCP4_DATA_PATH'
              value: '/mnt/ccp4data'
            }
            {
              name: 'CCP4'
              value: '/mnt/ccp4data/${ccp4Version}'
            }
            {
              name: 'LD_LIBRARY_PATH'
              value: '/mnt/ccp4data/${ccp4Version}/lib'
            }
            {
              name: 'CCP4I2_PROJECTS_DIR'
              value: '/mnt/projects'
            }
            {
              name: 'SERVICE_BUS_CONNECTION_STRING'
              secretRef: 'servicebus-connection'
            }
            {
              name: 'SERVICE_BUS_QUEUE_NAME'
              value: '${prefix}-jobs'
            }
            // Azure Storage for staged uploads (large file downloads in worker)
            {
              name: 'AZURE_STORAGE_ACCOUNT_NAME'
              value: storageAccountName
            }
            // Managed Identity Client ID for DefaultAzureCredential
            {
              name: 'AZURE_CLIENT_ID'
              value: containerAppsIdentityClientId
            }
          ]
          volumeMounts: concat(ccp4VolumeMount, [
            {
              volumeName: 'staticfiles-volume'
              mountPath: '/mnt/staticfiles'
            }
            {
              volumeName: 'mediafiles-volume'
              mountPath: '/mnt/mediafiles'
            }
            {
              volumeName: 'projects-volume'
              mountPath: '/mnt/projects'
            }
          ])
        }
      ]
      scale: {
        minReplicas: 0 // Scale to zero when no jobs - cost optimization for episodic workloads
        maxReplicas: 20  // Allow burst to 20 workers for parallel job processing
        rules: [
          {
            name: 'queue-scaling'
            custom: {
              type: 'azure-servicebus'
              metadata: {
                queueName: '${prefix}-jobs'
                namespace: '${prefix}-servicebus'
                messageCount: '5' // Scale up when 5+ messages in queue
              }
              auth: [
                {
                  secretRef: 'servicebus-connection'
                  triggerParameter: 'connection'
                }
              ]
            }
          }
        ]
      }
      volumes: concat(ccp4Volume, [
        {
          name: 'staticfiles-volume'
          storageName: 'staticfiles-private'
          storageType: 'AzureFile'
        }
        {
          name: 'mediafiles-volume'
          storageName: 'mediafiles-private'
          storageType: 'AzureFile'
        }
        {
          name: 'projects-volume'
          storageName: 'projects-private'
          storageType: 'AzureFile'
        }
      ])
    }
  }
}

// Web Container App
resource webApp 'Microsoft.App/containerApps@2023-05-01' = {
  name: webAppName
  location: resourceGroup().location
  properties: {
    managedEnvironmentId: containerAppsEnvironmentId
    configuration: {
      activeRevisionsMode: 'Single'
      ingress: {
        external: true
        targetPort: 3000
        allowInsecure: false
        traffic: [
          {
            weight: 100
            latestRevision: true
          }
        ]
      }
      registries: [
        {
          server: acrLoginServer
          username: acrName
          passwordSecretRef: 'registry-password'
        }
      ]
      secrets: [
        {
          name: 'registry-password'
          value: containerRegistry.listCredentials().passwords[0].value
        }
      ]
    }
    template: {
      containers: [
        {
          name: 'web'
          image: '${acrLoginServer}/ccp4i2/web:${imageTagWeb}'
          resources: {
            cpu: json('1.0')
            memory: '2.0Gi'
          }
          env: [
            {
              name: 'NEXT_PUBLIC_API_BASE_URL'
              value: 'http://${serverInternalFqdn}'
            }
            {
              name: 'API_BASE_URL'
              value: 'http://${serverInternalFqdn}'
            }
            {
              name: 'NEXT_PUBLIC_AAD_CLIENT_ID'
              value: aadClientId
            }
            {
              name: 'NEXT_PUBLIC_AAD_TENANT_ID'
              value: aadTenantId
            }
            {
              name: 'NEXT_PUBLIC_REQUIRE_AUTH'
              value: 'true'
            }
            // Enable proxy authentication (requires valid Azure AD token for API calls)
            {
              name: 'REQUIRE_PROXY_AUTH'
              value: 'true'
            }
          ]
          volumeMounts: concat(ccp4VolumeMount, [
            {
              volumeName: 'staticfiles-volume'
              mountPath: '/mnt/staticfiles'
            }
            {
              volumeName: 'mediafiles-volume'
              mountPath: '/mnt/mediafiles'
            }
            {
              volumeName: 'projects-volume'
              mountPath: '/mnt/projects'
            }
          ])
        }
      ]
      scale: {
        minReplicas: 1
        maxReplicas: 5
        rules: [
          {
            name: 'http-scaling'
            http: {
              metadata: {
                concurrentRequests: '50'
              }
            }
          }
        ]
      }
      volumes: concat(ccp4Volume, [
        {
          name: 'staticfiles-volume'
          storageName: 'staticfiles-private'
          storageType: 'AzureFile'
        }
        {
          name: 'mediafiles-volume'
          storageName: 'mediafiles-private'
          storageType: 'AzureFile'
        }
        {
          name: 'projects-volume'
          storageName: 'projects-private'
          storageType: 'AzureFile'
        }
      ])
    }
  }
}

// Authentication configuration for Server App removed - authentication now handled in frontend

// Authentication configuration for Web App removed - authentication now handled in frontend

// Key Vault RBAC Role Assignment for Server App (Key Vault Secrets User)
// NOTE: Role assignments are handled separately to avoid conflicts on redeployment
// resource keyVaultRoleAssignment 'Microsoft.Authorization/roleAssignments@2022-04-01' = {
//   name: guid(keyVault.id, serverApp.id, '4633458b-17de-408a-b874-0445c86b69e6', roleAssignmentSuffix)
//   scope: keyVault
//   properties: {
//     roleDefinitionId: subscriptionResourceId(
//       'Microsoft.Authorization/roleDefinitions',
//       '4633458b-17de-408a-b874-0445c86b69e6'
//     ) // Key Vault Secrets User
//     principalId: serverApp.identity.principalId
//     principalType: 'ServicePrincipal'
//   }
// }

// Key Vault RBAC Role Assignment for Worker App (Key Vault Secrets User)
// NOTE: Role assignments are handled separately to avoid conflicts on redeployment
// resource keyVaultRoleAssignmentWorker 'Microsoft.Authorization/roleAssignments@2022-04-01' = {
//   name: guid(keyVault.id, workerApp.id, '4633458b-17de-408a-b874-0445c86b69e6', roleAssignmentSuffix)
//   scope: keyVault
//   properties: {
//     roleDefinitionId: subscriptionResourceId(
//       'Microsoft.Authorization/roleDefinitions',
//       '4633458b-17de-408a-b874-0445c86b69e6'
//     ) // Key Vault Secrets User
//     principalId: workerApp.identity.principalId
//     principalType: 'ServicePrincipal'
//   }
// }

// Outputs
output serverInternalUrl string = 'http://${serverAppName}' // Internal only - no public URL
output webUrl string = 'https://${webApp.properties.configuration.ingress.fqdn}'
output workerAppName string = workerApp.name
