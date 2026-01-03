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

@description('Server container image tag')
param imageTagServer string = 'latest'

@description('Resource naming prefix')
param prefix string = 'ccp4i2-bicep'

@description('Shared Container Apps Identity ID')
param containerAppsIdentityId string

@description('Shared Container Apps Identity Principal ID')
param containerAppsIdentityPrincipalId string

// - PostgreSQL is accessed via private endpoint (no public access)
// - Key Vault is accessed via private endpoint (no public access)
// - Storage Account is accessed via private endpoint (no public access)
// - Container Apps run in dedicated VNet subnet with proper delegation
// - All communication happens within the Azure private network

// Management Container App Purpose:
// This is a long-lived container for interactive debugging and maintenance tasks.
// It runs 'tail -f /dev/null' to stay alive, allowing you to exec into it.
// Access via: az containerapp exec -g <rg> -n <app-name> --command bash
// Use cases:
//   - Interactive Django shell (python manage.py shell)
//   - Database migrations and inspection
//   - File system debugging (mounted volumes)
//   - Environment variable inspection

// Variables
var managementAppName = '${prefix}-management'
var webAppName = '${prefix}-web'

// Existing resources
resource containerRegistry 'Microsoft.ContainerRegistry/registries@2023-07-01' existing = {
  name: acrName
}

resource keyVault 'Microsoft.KeyVault/vaults@2023-02-01' existing = {
  name: keyVaultName
}

// Management Container App
resource managementApp 'Microsoft.App/containerApps@2023-05-01' = {
  name: managementAppName
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
      // No ingress - this is a management/debug container accessed via 'az containerapp exec'
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
          // Override default command to keep container running for interactive access
          command: ['/bin/bash']
          args: ['-c', 'export PYTHONPATH="/mnt/ccp4data/py-packages:$PYTHONPATH" && echo "Management container ready for interactive access" && tail -f /dev/null']
          resources: {
            cpu: json('2.0')
            memory: '4.0Gi'
          }
          // Remove health probes since we're not running the Django server
          // This is a management/debug container, not a production service
          env: [            {
              name: 'EXECUTION_MODE'
              value: 'azure'
            }
            {
              name: 'DJANGO_SETTINGS_MODULE'
              value: 'ccp4x.config.settings'
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
              value: 'true' // Ensure DEBUG is false in production
            }
            {
              name: 'CCP4_DATA_PATH'
              value: '/mnt/ccp4data'
            }
            {
              name: 'CCP4I2_PROJECTS_DIR'
              value: '/mnt/ccp4data/ccp4i2-projects'
            }
            {
              name: 'ALLOWED_HOSTS'
              value: '${managementAppName},${managementAppName}.whitecliff-258bc831.northeurope.azurecontainerapps.io,localhost,127.0.0.1,*'
            }
            {
              name: 'CORS_ALLOWED_ORIGINS'
              value: 'http://${webAppName},https://${webAppName}.whitecliff-258bc831.northeurope.azurecontainerapps.io'
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
          ]
          volumeMounts: [
            {
              volumeName: 'ccp4data-volume'
              mountPath: '/mnt/ccp4data'
            }
            {
              volumeName: 'staticfiles-volume'
              mountPath: '/mnt/staticfiles'
            }
            {
              volumeName: 'mediafiles-volume'
              mountPath: '/mnt/mediafiles'
            }
          ]
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
      volumes: [
        {
          name: 'ccp4data-volume'
          storageName: 'ccp4data-mount'
          storageType: 'AzureFile'
        }
        {
          name: 'staticfiles-volume'
          storageName: 'staticfiles-mount'
          storageType: 'AzureFile'
        }
        {
          name: 'mediafiles-volume'
          storageName: 'mediafiles-mount'
          storageType: 'AzureFile'
        }
      ]
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
output managementAppName string = managementApp.name
