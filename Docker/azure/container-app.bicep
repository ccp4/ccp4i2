// CCP4i2 Azure Container Apps Deployment
//
// This Bicep template deploys CCP4i2 to Azure Container Apps with:
// - Django API container
// - Next.js web client container
// - Azure PostgreSQL Flexible Server
// - Azure Files for CCP4 and projects storage
// - Application Insights monitoring
//
// Deploy:
//   az deployment group create \
//     --resource-group ccp4i2-rg \
//     --template-file Docker/azure/container-app.bicep \
//     --parameters @Docker/azure/parameters.json

@description('Base name for all resources')
param baseName string = 'ccp4i2'

@description('Location for all resources')
param location string = resourceGroup().location

@description('Container image for the server')
param serverImage string = 'ghcr.io/ccp4/ccp4i2-server:latest'

@description('Container image for the client')
param clientImage string = 'ghcr.io/ccp4/ccp4i2-client:latest'

@description('PostgreSQL administrator password')
@secure()
param dbPassword string

@description('Django secret key')
@secure()
param secretKey string

@description('CCP4 version directory name (e.g., ccp4-9, ccp4-20251105)')
param ccp4Version string = 'ccp4-9'

// Log Analytics Workspace
resource logAnalytics 'Microsoft.OperationalInsights/workspaces@2022-10-01' = {
  name: '${baseName}-logs'
  location: location
  properties: {
    sku: {
      name: 'PerGB2018'
    }
    retentionInDays: 30
  }
}

// Application Insights
resource appInsights 'Microsoft.Insights/components@2020-02-02' = {
  name: '${baseName}-insights'
  location: location
  kind: 'web'
  properties: {
    Application_Type: 'web'
    WorkspaceResourceId: logAnalytics.id
  }
}

// Container Apps Environment
resource containerAppEnv 'Microsoft.App/managedEnvironments@2023-05-01' = {
  name: '${baseName}-env'
  location: location
  properties: {
    appLogsConfiguration: {
      destination: 'log-analytics'
      logAnalyticsConfiguration: {
        customerId: logAnalytics.properties.customerId
        sharedKey: logAnalytics.listKeys().primarySharedKey
      }
    }
  }
}

// Storage Account for Azure Files
resource storageAccount 'Microsoft.Storage/storageAccounts@2023-01-01' = {
  name: '${baseName}storage'
  location: location
  sku: {
    name: 'Standard_LRS'
  }
  kind: 'StorageV2'
  properties: {
    accessTier: 'Hot'
  }
}

// File Share for CCP4 data
resource ccp4FileShare 'Microsoft.Storage/storageAccounts/fileServices/shares@2023-01-01' = {
  name: '${storageAccount.name}/default/ccp4data'
  properties: {
    shareQuota: 100  // GB
  }
}

// File Share for projects
resource projectsFileShare 'Microsoft.Storage/storageAccounts/fileServices/shares@2023-01-01' = {
  name: '${storageAccount.name}/default/projects'
  properties: {
    shareQuota: 50  // GB
  }
}

// Storage mount for Container Apps
resource storageMount 'Microsoft.App/managedEnvironments/storages@2023-05-01' = {
  name: 'ccp4storage'
  parent: containerAppEnv
  properties: {
    azureFile: {
      accountName: storageAccount.name
      accountKey: storageAccount.listKeys().keys[0].value
      shareName: 'ccp4data'
      accessMode: 'ReadOnly'
    }
  }
}

resource projectsMount 'Microsoft.App/managedEnvironments/storages@2023-05-01' = {
  name: 'projectsstorage'
  parent: containerAppEnv
  properties: {
    azureFile: {
      accountName: storageAccount.name
      accountKey: storageAccount.listKeys().keys[0].value
      shareName: 'projects'
      accessMode: 'ReadWrite'
    }
  }
}

// PostgreSQL Flexible Server
resource postgresServer 'Microsoft.DBforPostgreSQL/flexibleServers@2023-03-01-preview' = {
  name: '${baseName}-db'
  location: location
  sku: {
    name: 'Standard_B1ms'
    tier: 'Burstable'
  }
  properties: {
    version: '15'
    administratorLogin: 'ccp4admin'
    administratorLoginPassword: dbPassword
    storage: {
      storageSizeGB: 32
    }
    backup: {
      backupRetentionDays: 7
      geoRedundantBackup: 'Disabled'
    }
  }
}

// PostgreSQL Database
resource postgresDb 'Microsoft.DBforPostgreSQL/flexibleServers/databases@2023-03-01-preview' = {
  name: 'ccp4i2'
  parent: postgresServer
  properties: {
    charset: 'UTF8'
    collation: 'en_US.utf8'
  }
}

// Allow Azure services to access PostgreSQL
resource postgresFirewall 'Microsoft.DBforPostgreSQL/flexibleServers/firewallRules@2023-03-01-preview' = {
  name: 'AllowAzureServices'
  parent: postgresServer
  properties: {
    startIpAddress: '0.0.0.0'
    endIpAddress: '0.0.0.0'
  }
}

// Server Container App
resource serverApp 'Microsoft.App/containerApps@2023-05-01' = {
  name: '${baseName}-server'
  location: location
  properties: {
    managedEnvironmentId: containerAppEnv.id
    configuration: {
      ingress: {
        external: true
        targetPort: 8000
        transport: 'http'
        corsPolicy: {
          allowedOrigins: [
            'https://${baseName}-client.${containerAppEnv.properties.defaultDomain}'
          ]
          allowedMethods: ['GET', 'POST', 'PUT', 'DELETE', 'OPTIONS']
          allowedHeaders: ['*']
          allowCredentials: true
        }
      }
      secrets: [
        {
          name: 'db-password'
          value: dbPassword
        }
        {
          name: 'secret-key'
          value: secretKey
        }
      ]
    }
    template: {
      containers: [
        {
          name: 'server'
          image: serverImage
          resources: {
            cpu: json('0.5')
            memory: '1Gi'
          }
          env: [
            {
              name: 'SECRET_KEY'
              secretRef: 'secret-key'
            }
            {
              name: 'DB_HOST'
              value: postgresServer.properties.fullyQualifiedDomainName
            }
            {
              name: 'DB_USER'
              value: 'ccp4admin'
            }
            {
              name: 'DB_PASSWORD'
              secretRef: 'db-password'
            }
            {
              name: 'DB_NAME'
              value: 'ccp4i2'
            }
            {
              name: 'DJANGO_SETTINGS_MODULE'
              value: 'ccp4i2.config.settings'
            }
            {
              name: 'CCP4_DATA_PATH'
              value: '/mnt/ccp4data'
            }
            {
              name: 'CCP4_VERSION'
              value: ccp4Version
            }
            {
              name: 'CCP4I2_PROJECTS_DIR'
              value: '/mnt/projects'
            }
            {
              name: 'APPLICATIONINSIGHTS_CONNECTION_STRING'
              value: appInsights.properties.ConnectionString
            }
          ]
          volumeMounts: [
            {
              volumeName: 'ccp4data'
              mountPath: '/mnt/ccp4data'
            }
            {
              volumeName: 'projects'
              mountPath: '/mnt/projects'
            }
          ]
        }
      ]
      volumes: [
        {
          name: 'ccp4data'
          storageName: 'ccp4storage'
          storageType: 'AzureFile'
        }
        {
          name: 'projects'
          storageName: 'projectsstorage'
          storageType: 'AzureFile'
        }
      ]
      scale: {
        minReplicas: 1
        maxReplicas: 3
        rules: [
          {
            name: 'http-rule'
            http: {
              metadata: {
                concurrentRequests: '50'
              }
            }
          }
        ]
      }
    }
  }
}

// Client Container App
resource clientApp 'Microsoft.App/containerApps@2023-05-01' = {
  name: '${baseName}-client'
  location: location
  properties: {
    managedEnvironmentId: containerAppEnv.id
    configuration: {
      ingress: {
        external: true
        targetPort: 3000
        transport: 'http'
      }
    }
    template: {
      containers: [
        {
          name: 'client'
          image: clientImage
          resources: {
            cpu: json('0.25')
            memory: '0.5Gi'
          }
          env: [
            {
              name: 'NEXT_PUBLIC_API_URL'
              value: 'https://${serverApp.properties.configuration.ingress.fqdn}'
            }
          ]
        }
      ]
      scale: {
        minReplicas: 1
        maxReplicas: 3
      }
    }
  }
}

// Outputs
output serverUrl string = 'https://${serverApp.properties.configuration.ingress.fqdn}'
output clientUrl string = 'https://${clientApp.properties.configuration.ingress.fqdn}'
output storageAccountName string = storageAccount.name
output postgresServerName string = postgresServer.name
