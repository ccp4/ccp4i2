// Self-contained Ubuntu 24.04 desktop test VM for exercising the CCP4i2 Electron
// app (AppImage + .deb) on a MODERN, locked-down kernel.
//
// Why 24.04: it ships kernel.apparmor_restrict_unprivileged_userns=1 by default
// — the exact restriction that breaks the AppImage's Chromium sandbox and that
// 22.04 simply does not have. A 22.04 rig green-lights everything and reproduces
// nothing.
//
// Deliberately self-contained (own VNet/subnet, no shared/private network): this
// is a throwaway rig, so the whole resource group can be created and deleted
// freely without touching any other infrastructure.

@description('Location for all resources')
param location string = resourceGroup().location

@description('Resource naming prefix')
param prefix string = 'ccp4i2-linux-test'

@description('VM administrator username')
param adminUsername string = 'ccp4admin'

@description('VM administrator password (used for RDP/xrdp and password SSH)')
@secure()
param adminPassword string

@description('SSH public key for key-based SSH login')
param adminSshPublicKey string

@description('VM size - D4s_v5 gives 4 vCPU, 16GB RAM')
param vmSize string = 'Standard_D4s_v5'

@description('Auto-shutdown time in GMT (HHMM, e.g. 1900 = 7pm) to avoid burning money')
param autoShutdownTime string = '1900'

@description('Your IP address (CIDR) for SSH + RDP access, e.g. 1.2.3.4/32')
param allowedSourceIp string

var vmName = '${prefix}-vm'
var nicName = '${vmName}-nic'
var publicIpName = '${vmName}-pip'
var nsgName = '${vmName}-nsg'
var osDiskName = '${vmName}-osdisk'
var vnetName = '${prefix}-vnet'
var subnetName = 'default'

// SSH + RDP restricted to your IP
resource nsg 'Microsoft.Network/networkSecurityGroups@2023-05-01' = {
  name: nsgName
  location: location
  properties: {
    securityRules: [
      {
        name: 'AllowSSH'
        properties: {
          protocol: 'Tcp'
          sourcePortRange: '*'
          destinationPortRange: '22'
          sourceAddressPrefix: allowedSourceIp
          destinationAddressPrefix: '*'
          access: 'Allow'
          priority: 1000
          direction: 'Inbound'
        }
      }
      {
        name: 'AllowRDP'
        properties: {
          protocol: 'Tcp'
          sourcePortRange: '*'
          destinationPortRange: '3389'
          sourceAddressPrefix: allowedSourceIp
          destinationAddressPrefix: '*'
          access: 'Allow'
          priority: 1010
          direction: 'Inbound'
        }
      }
    ]
  }
}

// Own VNet + subnet — no dependency on any shared network.
resource vnet 'Microsoft.Network/virtualNetworks@2023-05-01' = {
  name: vnetName
  location: location
  properties: {
    addressSpace: {
      addressPrefixes: [ '10.42.0.0/24' ]
    }
    subnets: [
      {
        name: subnetName
        properties: {
          addressPrefix: '10.42.0.0/24'
          networkSecurityGroup: {
            id: nsg.id
          }
        }
      }
    ]
  }
}

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
            id: resourceId('Microsoft.Network/virtualNetworks/subnets', vnetName, subnetName)
          }
        }
      }
    ]
  }
  dependsOn: [
    vnet
  ]
}

// Ubuntu 24.04 LTS (Noble) — the modern locked-down kernel we need to reproduce.
resource vm 'Microsoft.Compute/virtualMachines@2023-07-01' = {
  name: vmName
  location: location
  properties: {
    hardwareProfile: {
      vmSize: vmSize
    }
    osProfile: {
      computerName: 'ccp4i2linux'
      adminUsername: adminUsername
      adminPassword: adminPassword
      // Password auth stays ON (xrdp login needs it); SSH key added too.
      linuxConfiguration: {
        disablePasswordAuthentication: false
        ssh: {
          publicKeys: [
            {
              path: '/home/${adminUsername}/.ssh/authorized_keys'
              keyData: adminSshPublicKey
            }
          ]
        }
      }
      customData: loadFileAsBase64('linux-test-vm.cloud-init.yaml')
    }
    storageProfile: {
      imageReference: {
        publisher: 'Canonical'
        offer: 'ubuntu-24_04-lts'
        // For the 24.04 offer, 'server' is the Gen2 image ('server-gen1' is Gen1);
        // the old '<ver>-lts-gen2' convention doesn't apply here.
        sku: 'server'
        version: 'latest'
      }
      osDisk: {
        name: osDiskName
        createOption: 'FromImage'
        managedDisk: {
          storageAccountType: 'Premium_LRS'
        }
        // CCP4 extracted (~15GB) + tarball (~6GB) + headroom
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

// Auto-shutdown to avoid burning money if left running.
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

output vmName string = vm.name
output publicIpAddress string = publicIp.properties.ipAddress
output sshCommand string = 'ssh ${adminUsername}@${publicIp.properties.ipAddress}'
output adminUsername string = adminUsername
