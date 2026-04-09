/**
 * Teams API Diagnostic Component
 * Temporary component to test Microsoft Teams integration independently
 */

import React, { useState } from "react";
import { Button, Box, Typography, Paper, Alert, Divider } from "@mui/material";
import { useMsal } from "@azure/msal-react";
import {
  checkTeamsMembership,
  DEFAULT_TEAMS_CONFIG,
  MsalAuthProvider,
} from "../utils/teams-auth";
import { Client } from "@microsoft/microsoft-graph-client";

interface TestResult {
  step: string;
  status: "success" | "error" | "info";
  message: string;
  data?: any;
}

export const TeamsDiagnostic: React.FC = () => {
  const { instance, accounts } = useMsal();
  const [testing, setTesting] = useState(false);
  const [results, setResults] = useState<TestResult[]>([]);

  const addResult = (result: TestResult) => {
    setResults((prev) => [...prev, result]);
  };

  const clearResults = () => {
    setResults([]);
  };

  const requestTeamsPermission = async () => {
    try {
      addResult({
        step: "Teams Permission Request",
        status: "info",
        message: "Attempting to request Teams permission via redirect...",
      });

      // This will redirect to Microsoft consent page
      await instance.acquireTokenRedirect({
        scopes: ["User.Read", "Team.ReadBasic.All"],
        prompt: "consent",
      });
    } catch (error) {
      addResult({
        step: "Teams Permission Request",
        status: "error",
        message: `Failed to request Teams permission: ${error.message}`,
      });
    }
  };

  const generateAdminConsentInfo = () => {
    if (accounts.length === 0) return null;

    const account = accounts[0];

    // Get client ID from MSAL configuration
    const clientId = instance.getConfiguration().auth.clientId;
    const tenantId = account.tenantId;

    const adminConsentUrl = `https://login.microsoftonline.com/${tenantId}/adminconsent?client_id=${clientId}&redirect_uri=${encodeURIComponent(window.location.origin)}&scope=User.Read%20Team.ReadBasic.All`;

    addResult({
      step: "Admin Consent Required",
      status: "error",
      message:
        "Your organization requires admin approval for Teams permissions. Contact your IT administrator with the information below.",
      data: {
        adminConsentUrl,
        clientId,
        tenantId,
        requiredPermissions: ["User.Read", "Team.ReadBasic.All"],
        instructions: [
          "1. Send the admin consent URL to your IT administrator",
          "2. Ask them to grant 'Team.ReadBasic.All' permission for the CCP4i2 application",
          "3. Once approved, all Newcastle University users will have Teams authorization",
          "4. The permission only reads basic team membership - no sensitive data",
        ],
      },
    });
  };

  const testTeamsIntegration = async () => {
    setTesting(true);
    clearResults();

    try {
      // Step 1: Check if user is authenticated
      addResult({
        step: "Authentication Check",
        status: "info",
        message: `Found ${accounts.length} authenticated account(s)`,
        data: accounts.map((acc) => ({
          username: acc.username,
          tenantId: acc.tenantId,
        })),
      });

      if (accounts.length === 0) {
        addResult({
          step: "Authentication Check",
          status: "error",
          message: "No authenticated accounts found. Please log in first.",
        });
        return;
      }

      // Step 2: Test basic token acquisition
      try {
        const basicTokenResponse = await instance.acquireTokenSilent({
          scopes: ["User.Read"],
          account: accounts[0],
        });

        addResult({
          step: "Basic Token Acquisition",
          status: "success",
          message: "Successfully acquired User.Read token",
          data: {
            scopes: basicTokenResponse.scopes,
            expiresOn: basicTokenResponse.expiresOn,
          },
        });
      } catch (basicError) {
        addResult({
          step: "Basic Token Acquisition",
          status: "error",
          message: `Failed to acquire basic token: ${basicError.message}`,
        });
        return;
      }

      // Step 3: Test Teams token acquisition (silent)
      try {
        const teamsTokenResponse = await instance.acquireTokenSilent({
          scopes: ["User.Read", "Team.ReadBasic.All"],
          account: accounts[0],
        });

        addResult({
          step: "Teams Token (Silent)",
          status: "success",
          message: "Successfully acquired Teams token silently",
          data: {
            scopes: teamsTokenResponse.scopes,
            expiresOn: teamsTokenResponse.expiresOn,
          },
        });

        // Step 4: Test network connectivity first
        await testNetworkConnectivity();

        // Step 5: Test Graph API connection
        await testGraphApiConnection(teamsTokenResponse.accessToken);
      } catch (silentError) {
        addResult({
          step: "Teams Token (Silent)",
          status: "error",
          message: `Silent token acquisition failed: ${silentError.message}`,
          data: {
            errorCode: silentError.errorCode,
            errorMessage: silentError.errorMessage,
          },
        });

        // Step 3b: Try interactive Teams token (redirect instead of popup)
        try {
          addResult({
            step: "Teams Token (Interactive)",
            status: "info",
            message: "Popup failed, trying redirect-based consent...",
          });

          // Use redirect instead of popup for Electron compatibility
          await instance.acquireTokenRedirect({
            scopes: ["User.Read", "Team.ReadBasic.All"],
            prompt: "consent",
          });

          // Note: This will cause a page redirect, so we won't reach this point
          addResult({
            step: "Teams Token (Interactive)",
            status: "info",
            message:
              "Redirect initiated - you'll be taken to Microsoft consent page...",
          });
        } catch (redirectError) {
          addResult({
            step: "Teams Token (Interactive)",
            status: "error",
            message: `Redirect token acquisition failed: ${redirectError.message}`,
            data: {
              errorCode: redirectError.errorCode,
              errorMessage: redirectError.errorMessage,
            },
          });

          // Step 3c: Try alternative approach - check what scopes are actually available
          addResult({
            step: "Available Scopes Check",
            status: "info",
            message: "Checking what permissions are currently available...",
          });

          try {
            const accounts = instance.getAllAccounts();
            if (accounts.length > 0) {
              const account = accounts[0];
              addResult({
                step: "Account Details",
                status: "info",
                message: "Current account information",
                data: {
                  username: account.username,
                  tenantId: account.tenantId,
                  scopes: account.idTokenClaims?.scp || "No scopes in ID token",
                },
              });
            }
          } catch (accountError) {
            addResult({
              step: "Account Details",
              status: "error",
              message: `Failed to get account details: ${accountError.message}`,
            });
          }
        }
      }

      // Step 5: Test full Teams membership check
      try {
        addResult({
          step: "Teams Membership Check",
          status: "info",
          message: "Testing full Teams membership function...",
        });

        const teamsResult = await checkTeamsMembership(
          instance,
          DEFAULT_TEAMS_CONFIG
        );

        addResult({
          step: "Teams Membership Check",
          status: teamsResult.hasAccess ? "success" : "error",
          message: teamsResult.reason,
          data: {
            hasAccess: teamsResult.hasAccess,
            teamsFound: teamsResult.teamsFound,
            error: teamsResult.error,
          },
        });
      } catch (teamsError) {
        addResult({
          step: "Teams Membership Check",
          status: "error",
          message: `Teams membership check failed: ${teamsError.message}`,
        });
      }
    } catch (error) {
      addResult({
        step: "General Error",
        status: "error",
        message: `Unexpected error: ${error.message}`,
      });
    } finally {
      setTesting(false);
    }
  };

  const testNetworkConnectivity = async () => {
    try {
      addResult({
        step: "Network Connectivity",
        status: "info",
        message: "Testing basic network connectivity...",
      });

      // Test 1: Basic internet connectivity
      const googleResponse = await fetch("https://www.google.com", {
        mode: "no-cors",
        method: "HEAD",
      });

      addResult({
        step: "Internet Connectivity",
        status: "success",
        message: "Basic internet connectivity confirmed",
      });

      // Test 2: Microsoft Graph endpoint accessibility
      try {
        const graphResponse = await fetch("https://graph.microsoft.com/v1.0/", {
          mode: "no-cors",
          method: "HEAD",
        });

        addResult({
          step: "Graph Endpoint Accessibility",
          status: "success",
          message: "Microsoft Graph endpoint is accessible",
        });
      } catch (graphEndpointError) {
        addResult({
          step: "Graph Endpoint Accessibility",
          status: "error",
          message: `Graph endpoint test failed: ${graphEndpointError.message}`,
          data: {
            error: graphEndpointError.message,
            name: graphEndpointError.name,
          },
        });
      }
    } catch (networkError) {
      addResult({
        step: "Network Connectivity",
        status: "error",
        message: `Network connectivity test failed: ${networkError.message}`,
        data: {
          error: networkError.message,
          name: networkError.name,
        },
      });
    }
  };

  const testGraphApiConnection = async (accessToken: string) => {
    try {
      // Test 1: Direct fetch to Graph API (bypass Microsoft Graph SDK)
      try {
        addResult({
          step: "Direct Graph API Test",
          status: "info",
          message: "Testing direct fetch to Microsoft Graph API...",
        });

        const directResponse = await fetch(
          "https://graph.microsoft.com/v1.0/me",
          {
            method: "GET",
            headers: {
              Authorization: `Bearer ${accessToken}`,
              "Content-Type": "application/json",
            },
          }
        );

        if (!directResponse.ok) {
          throw new Error(
            `HTTP ${directResponse.status}: ${directResponse.statusText}`
          );
        }

        const userInfo = await directResponse.json();
        addResult({
          step: "Direct Graph API - User Info",
          status: "success",
          message: `Successfully retrieved user info for ${userInfo.displayName}`,
          data: {
            displayName: userInfo.displayName,
            mail: userInfo.mail,
            id: userInfo.id,
          },
        });

        // Test 2: Direct fetch for teams
        try {
          const teamsResponse = await fetch(
            "https://graph.microsoft.com/v1.0/me/joinedTeams",
            {
              method: "GET",
              headers: {
                Authorization: `Bearer ${accessToken}`,
                "Content-Type": "application/json",
              },
            }
          );

          if (!teamsResponse.ok) {
            throw new Error(
              `HTTP ${teamsResponse.status}: ${teamsResponse.statusText}`
            );
          }

          const teamsData = await teamsResponse.json();
          addResult({
            step: "Direct Graph API - User Teams",
            status: "success",
            message: `Found ${teamsData.value?.length || 0} team(s)`,
            data: teamsData.value?.map((team: any) => ({
              id: team.id,
              displayName: team.displayName,
              description: team.description,
            })),
          });
        } catch (teamsDirectError) {
          addResult({
            step: "Direct Graph API - User Teams",
            status: "error",
            message: `Direct teams API failed: ${teamsDirectError.message}`,
            data: {
              error: teamsDirectError.message,
              name: teamsDirectError.name,
              stack: teamsDirectError.stack?.substring(0, 500),
            },
          });
        }
      } catch (directError) {
        addResult({
          step: "Direct Graph API Test",
          status: "error",
          message: `Direct Graph API failed: ${directError.message}`,
          data: {
            error: directError.message,
            name: directError.name,
            stack: directError.stack?.substring(0, 500),
          },
        });

        // Fallback: Try with Microsoft Graph SDK
        addResult({
          step: "Microsoft Graph SDK Test",
          status: "info",
          message: "Direct API failed, trying Microsoft Graph SDK...",
        });

        try {
          const authProvider = new MsalAuthProvider(instance);
          const graphClient = Client.initWithMiddleware({ authProvider });

          const userInfo = await graphClient.api("/me").get();
          addResult({
            step: "Graph SDK - User Info",
            status: "success",
            message: `SDK successfully retrieved user info for ${userInfo.displayName}`,
            data: {
              displayName: userInfo.displayName,
              mail: userInfo.mail,
              id: userInfo.id,
            },
          });
        } catch (sdkError) {
          addResult({
            step: "Graph SDK Test",
            status: "error",
            message: `Microsoft Graph SDK failed: ${sdkError.message}`,
            data: {
              error: sdkError.message,
              name: sdkError.name,
              code: sdkError.code,
            },
          });
        }
      }
    } catch (graphError) {
      addResult({
        step: "Graph API Connection",
        status: "error",
        message: `Graph API connection failed: ${graphError.message}`,
        data: {
          error: graphError.message,
          name: graphError.name,
          stack: graphError.stack?.substring(0, 500),
        },
      });
    }
  };

  const getStatusColor = (status: string) => {
    switch (status) {
      case "success":
        return "success";
      case "error":
        return "error";
      case "info":
        return "info";
      default:
        return "info";
    }
  };

  return (
    <Box>
      <Typography variant="h6" gutterBottom>
        Teams API Diagnostic Tool
      </Typography>

      <Typography variant="body2" color="textSecondary" paragraph>
        This tool tests Microsoft Teams integration step by step to diagnose any
        issues.
      </Typography>

      <Alert severity="warning" sx={{ mb: 2 }}>
        <Typography variant="body2">
          <strong>Admin Consent Required:</strong> Based on your testing,
          Newcastle University requires admin approval for Teams permissions.
          Use "Generate Admin Consent Info" to get the URL and instructions for
          your IT administrator.
        </Typography>
      </Alert>

      <Alert severity="info" sx={{ mb: 2 }}>
        <Typography variant="body2">
          <strong>Popup Issues:</strong> If you see "popup_window_error", this
          is normal in Electron apps. The redirect method or admin consent
          approach will work better.
        </Typography>
      </Alert>

      <Box sx={{ mb: 2 }}>
        <Button
          variant="contained"
          onClick={testTeamsIntegration}
          disabled={testing}
          sx={{ mr: 2 }}
        >
          {testing ? "Running Tests..." : "Run Teams Diagnostic"}
        </Button>

        <Button
          variant="contained"
          color="secondary"
          onClick={requestTeamsPermission}
          disabled={testing}
          sx={{ mr: 2 }}
        >
          Request Teams Permission
        </Button>

        <Button
          variant="contained"
          color="warning"
          onClick={generateAdminConsentInfo}
          disabled={testing}
          sx={{ mr: 2 }}
        >
          Generate Admin Consent Info
        </Button>

        <Button variant="outlined" onClick={clearResults} disabled={testing}>
          Clear Results
        </Button>
      </Box>

      {results.length > 0 && (
        <Box sx={{ mt: 3 }}>
          <Typography variant="h6" gutterBottom>
            Test Results:
          </Typography>

          {results.map((result, index) => (
            <Box key={index} sx={{ mb: 2 }}>
              <Alert severity={getStatusColor(result.status)}>
                <Typography variant="subtitle2" component="div">
                  {result.step}
                </Typography>
                <Typography variant="body2">{result.message}</Typography>
                {result.data && (
                  <Box sx={{ mt: 1 }}>
                    <Typography variant="caption" component="div">
                      Data:
                    </Typography>
                    <pre
                      style={{
                        fontSize: "0.75rem",
                        margin: 0,
                        whiteSpace: "pre-wrap",
                      }}
                    >
                      {JSON.stringify(result.data, null, 2)}
                    </pre>
                  </Box>
                )}
              </Alert>
            </Box>
          ))}
        </Box>
      )}
    </Box>
  );
};
