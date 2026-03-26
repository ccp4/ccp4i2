@description('The location for all resources')
param location string = resourceGroup().location

@description('Resource naming prefix')
param prefix string = 'ccp4i2-bicep'

@description('Environment suffix')
param environment string = 'uk'

@description('VM administrator username')
param adminUsername string = 'ccp4admin'

@description('VM administrator password')
@secure()
param adminPassword string

@description('VM size - D4s_v5 gives 4 vCPU, 16GB RAM')
param vmSize string = 'Standard_D4s_v5'

@description('Auto-shutdown time in UTC (HHMM format, e.g. 1900 = 7pm)')
param autoShutdownTime string = '1900'

@description('Your IP address for RDP access (use /32 for single IP)')
param allowedRdpSourceIp string

// Variables
var vmName = '${prefix}-win-test'
var nicName = '${vmName}-nic'
var publicIpName = '${vmName}-pip'
var nsgName = '${vmName}-nsg'
var osDiskName = '${vmName}-osdisk'

// Network Security Group - RDP access restricted to your IP
resource nsg 'Microsoft.Network/networkSecurityGroups@2023-05-01' = {
  name: nsgName
  location: location
  properties: {
    securityRules: [
      {
        name: 'AllowRDP'
        properties: {
          protocol: 'Tcp'
          sourcePortRange: '*'
          destinationPortRange: '3389'
          sourceAddressPrefix: allowedRdpSourceIp
          destinationAddressPrefix: '*'
          access: 'Allow'
          priority: 1000
          direction: 'Inbound'
        }
      }
    ]
  }
}

// Public IP for RDP access
resource publicIp 'Microsoft.Network/publicIPAddresses@2023-05-01' = {
  name: publicIpName
  location: location
  sku: {
    name: 'Standard'
  }
  properties: {
    publicIPAllocationMethod: 'Static'
  }
}

// Reference existing VNet from infrastructure.bicep
resource existingVnet 'Microsoft.Network/virtualNetworks@2023-05-01' existing = {
  name: '${prefix}-vnet-${environment}'
}

// Add a subnet for the test VM (10.0.4.0/24 - next available range)
resource vmSubnet 'Microsoft.Network/virtualNetworks/subnets@2023-05-01' = {
  parent: existingVnet
  name: 'test-vm-subnet'
  properties: {
    addressPrefix: '10.0.4.0/24'
    networkSecurityGroup: {
      id: nsg.id
    }
  }
}

// Network Interface
resource nic 'Microsoft.Network/networkInterfaces@2023-05-01' = {
  name: nicName
  location: location
  properties: {
    ipConfigurations: [
      {
        name: 'ipconfig1'
        properties: {
          privateIPAllocationMethod: 'Dynamic'
          publicIPAddress: {
            id: publicIp.id
          }
          subnet: {
            id: vmSubnet.id
          }
        }
      }
    ]
  }
}

// Windows Server 2022 VM
resource vm 'Microsoft.Compute/virtualMachines@2023-07-01' = {
  name: vmName
  location: location
  properties: {
    hardwareProfile: {
      vmSize: vmSize
    }
    osProfile: {
      computerName: 'ccp4i2test'
      adminUsername: adminUsername
      adminPassword: adminPassword
    }
    storageProfile: {
      imageReference: {
        publisher: 'MicrosoftWindowsServer'
        offer: 'WindowsServer'
        sku: '2022-datacenter-g2'
        version: 'latest'
      }
      osDisk: {
        name: osDiskName
        createOption: 'FromImage'
        managedDisk: {
          storageAccountType: 'Premium_LRS'
        }
        diskSizeGB: 256
      }
    }
    networkProfile: {
      networkInterfaces: [
        {
          id: nic.id
        }
      ]
    }
  }
}

// Auto-shutdown schedule to avoid burning money
resource autoShutdown 'Microsoft.DevTestLab/schedules@2018-09-15' = {
  name: 'shutdown-computevm-${vmName}'
  location: location
  properties: {
    status: 'Enabled'
    taskType: 'ComputeVmShutdownTask'
    dailyRecurrence: {
      time: autoShutdownTime
    }
    timeZoneId: 'GMT Standard Time'
    targetResourceId: vm.id
    notificationSettings: {
      status: 'Disabled'
    }
  }
}

// Outputs
output vmName string = vm.name
output publicIpAddress string = publicIp.properties.ipAddress
output rdpCommand string = 'mstsc /v:${publicIp.properties.ipAddress}'
output adminUsername string = adminUsername
