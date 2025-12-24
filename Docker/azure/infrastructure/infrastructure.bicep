@description('The location for all resources')
param location string = resourceGroup().location

@description('Resource naming prefix')
param prefix string = 'ccp4i2-bicep'

@description('Environment suffix (e.g., dev, staging, prod)')
param environment string = 'ne'

@description('PostgreSQL administrator password')
@secure()
param postgresAdminPassword string

// Variables
var resourceSuffix = '${prefix}-${environment}'
var acrName = 'ccp4acr${environment}${substring(uniqueString(resourceGroup().id), 0, 4)}' // Keep it short and unique
var storageAccountName = 'stor${environment}${substring(uniqueString(resourceGroup().id), 0, 8)}' // Shorter for 24 char limit
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

// Storage Account
resource storageAccount 'Microsoft.Storage/storageAccounts@2023-01-01' = {
  name: storageAccountName
  location: location
  sku: {
    name: 'Standard_LRS'
  }
  kind: 'StorageV2'
  properties: {
    allowBlobPublicAccess: false
    allowSharedKeyAccess: true
    minimumTlsVersion: 'TLS1_2'
    supportsHttpsTrafficOnly: true
    publicNetworkAccess: 'Disabled'
    networkAcls: {
      defaultAction: 'Deny'
      bypass: 'AzureServices'
    }
  }
}

// File Shares
resource ccp4dataShare 'Microsoft.Storage/storageAccounts/fileServices/shares@2023-01-01' = {
  name: '${storageAccount.name}/default/ccp4data'
  properties: {
    shareQuota: 100
  }
}

resource staticfilesShare 'Microsoft.Storage/storageAccounts/fileServices/shares@2023-01-01' = {
  name: '${storageAccount.name}/default/staticfiles'
  properties: {
    shareQuota: 10
  }
}

resource mediafilesShare 'Microsoft.Storage/storageAccounts/fileServices/shares@2023-01-01' = {
  name: '${storageAccount.name}/default/mediafiles'
  properties: {
    shareQuota: 50
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

// PostgreSQL Flexible Server
resource postgresServer 'Microsoft.DBforPostgreSQL/flexibleServers@2023-03-01-preview' = {
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

// Remove the old firewall rule since we're using private endpoints

// Private Endpoints
resource storagePrivateEndpoint 'Microsoft.Network/privateEndpoints@2023-05-01' = {
  name: '${resourceSuffix}-storage-pe'
  location: location
  properties: {
    subnet: {
      id: '${vnet.id}/subnets/private-endpoints-subnet'
    }
    privateLinkServiceConnections: [
      {
        name: 'storage-connection'
        properties: {
          privateLinkServiceId: storageAccount.id
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
          privateLinkServiceId: postgresServer.id
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

// Storage for Container Apps Environment
resource containerAppsStorage 'Microsoft.App/managedEnvironments/storages@2023-05-01' = {
  name: 'ccp4data-mount'
  parent: containerAppsEnvironment
  properties: {
    azureFile: {
      accountName: storageAccount.name
      accountKey: storageAccount.listKeys().keys[0].value
      shareName: 'ccp4data'
      accessMode: 'ReadWrite'
    }
  }
}

resource containerAppsStaticStorage 'Microsoft.App/managedEnvironments/storages@2023-05-01' = {
  name: 'staticfiles-mount'
  parent: containerAppsEnvironment
  properties: {
    azureFile: {
      accountName: storageAccount.name
      accountKey: storageAccount.listKeys().keys[0].value
      shareName: 'staticfiles'
      accessMode: 'ReadWrite'
    }
  }
}

resource containerAppsMediaStorage 'Microsoft.App/managedEnvironments/storages@2023-05-01' = {
  name: 'mediafiles-mount'
  parent: containerAppsEnvironment
  properties: {
    azureFile: {
      accountName: storageAccount.name
      accountKey: storageAccount.listKeys().keys[0].value
      shareName: 'mediafiles'
      accessMode: 'ReadWrite'
    }
  }
}

// Outputs
output acrName string = containerRegistry.name
output acrLoginServer string = containerRegistry.properties.loginServer
output storageAccountName string = storageAccount.name
output keyVaultName string = keyVault.name
output postgresServerName string = postgresServer.name
output postgresServerFqdn string = postgresServer.properties.fullyQualifiedDomainName
output containerAppsEnvironmentName string = containerAppsEnvironment.name
output containerAppsEnvironmentId string = containerAppsEnvironment.id
output vnetId string = vnet.id
output vnetName string = vnet.name
output containerAppsSubnetId string = '${vnet.id}/subnets/container-apps-subnet'
output privateEndpointsSubnetId string = '${vnet.id}/subnets/private-endpoints-subnet'
output resourceGroupName string = resourceGroup().name
output containerAppsIdentityId string = containerAppsIdentity.id
output containerAppsIdentityPrincipalId string = containerAppsIdentity.properties.principalId
