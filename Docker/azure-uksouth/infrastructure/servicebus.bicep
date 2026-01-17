@description('Resource group location')
param location string = resourceGroup().location

@description('Prefix for resource naming')
param prefix string = 'ccp4i2-bicep'

@description('Key Vault name')
param keyVaultName string

@description('Service Bus SKU')
param sbSku string = 'Premium'

var sbNamespaceName = '${prefix}-servicebus'
var sbQueueName = '${prefix}-jobs'

resource keyVault 'Microsoft.KeyVault/vaults@2023-02-01' existing = {
  name: keyVaultName
}

resource sbNamespace 'Microsoft.ServiceBus/namespaces@2024-01-01' = {
  name: sbNamespaceName
  location: location
  sku: {
    name: sbSku
    tier: sbSku
  }
  properties: {
    zoneRedundant: false
    minimumTlsVersion: '1.2'
    publicNetworkAccess: 'Disabled' // Use private endpoints for security
  }
}

resource sbQueue 'Microsoft.ServiceBus/namespaces/queues@2024-01-01' = {
  name: '${sbNamespaceName}/${sbQueueName}'
  properties: {
    enablePartitioning: true
    maxSizeInMegabytes: 1024
    requiresDuplicateDetection: false
    enableBatchedOperations: true
    deadLetteringOnMessageExpiration: true
  }
}

resource sbAuthRule 'Microsoft.ServiceBus/namespaces/queues/authorizationRules@2024-01-01' = {
  name: '${sbNamespaceName}/${sbQueueName}/SendListen'
  properties: {
    rights: [
      'Send'
      'Listen'
    ]
  }
}

resource sbSecret 'Microsoft.KeyVault/vaults/secrets@2023-02-01' = {
  name: '${keyVault.name}/servicebus-connection'
  properties: {
    value: listKeys(sbAuthRule.id, sbAuthRule.apiVersion).primaryConnectionString
  }
}

output serviceBusNamespace string = sbNamespace.name
output serviceBusQueue string = sbQueue.name
output serviceBusConnectionSecret string = sbSecret.name
