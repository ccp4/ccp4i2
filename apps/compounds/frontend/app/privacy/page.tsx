'use client';

import { Container, Typography, Box, Paper, Link as MuiLink } from '@mui/material';
import Link from 'next/link';
import ArrowBackIcon from '@mui/icons-material/ArrowBack';

export default function PrivacyPolicyPage() {
  return (
    <Container maxWidth="md" sx={{ py: 4 }}>
      <Box sx={{ mb: 3 }}>
        <MuiLink
          component={Link}
          href="/"
          sx={{ display: 'inline-flex', alignItems: 'center', gap: 0.5 }}
        >
          <ArrowBackIcon fontSize="small" />
          Back to Home
        </MuiLink>
      </Box>

      <Paper elevation={1} sx={{ p: 4 }}>
        <Typography variant="h3" component="h1" gutterBottom>
          Privacy Policy
        </Typography>

        <Typography variant="body2" color="text.secondary" sx={{ mb: 4 }}>
          Last updated: January 2025
        </Typography>

        <Box sx={{ '& h2': { mt: 4, mb: 2 }, '& p': { mb: 2 } }}>
          <Typography variant="h5" component="h2">
            Introduction
          </Typography>
          <Typography>
            CCP4i2 is a crystallographic computing platform developed and operated by
            Newcastle University. This privacy policy explains how we collect, use, and
            protect your information when you use this service.
          </Typography>

          <Typography variant="h5" component="h2">
            Information We Collect
          </Typography>
          <Typography>
            <strong>Authentication Data:</strong> We use Microsoft Entra ID (Azure Active Directory)
            for authentication. When you sign in, we receive your name, email address, and
            organisational affiliation from your identity provider.
          </Typography>
          <Typography>
            <strong>Scientific Data:</strong> You may upload crystallographic data, molecular
            structures, compound information, and related scientific data. This data is stored
            to provide the service functionality.
          </Typography>
          <Typography>
            <strong>Usage Data:</strong> We collect logs of service usage for operational purposes,
            including job execution records and access timestamps.
          </Typography>

          <Typography variant="h5" component="h2">
            How We Use Your Information
          </Typography>
          <Typography>
            Your information is used to:
          </Typography>
          <Box component="ul" sx={{ pl: 3, mb: 2 }}>
            <li>Provide access to crystallographic computing tools and pipelines</li>
            <li>Store and manage your projects and scientific data</li>
            <li>Track job execution and provide results</li>
            <li>Maintain service security and diagnose technical issues</li>
          </Box>

          <Typography variant="h5" component="h2">
            Data Storage and Security
          </Typography>
          <Typography>
            Your data is stored on Microsoft Azure infrastructure located in the UK
            (UK South region). We implement appropriate technical and organisational
            measures to protect your data, including encryption in transit and at rest.
          </Typography>

          <Typography variant="h5" component="h2">
            Data Sharing
          </Typography>
          <Typography>
            We do not sell or share your personal or scientific data with third parties
            for commercial purposes. Data may be shared only:
          </Typography>
          <Box component="ul" sx={{ pl: 3, mb: 2 }}>
            <li>With your explicit consent</li>
            <li>To comply with legal obligations</li>
            <li>With service providers necessary for platform operation (Microsoft Azure)</li>
          </Box>

          <Typography variant="h5" component="h2">
            Data Retention
          </Typography>
          <Typography>
            Project data and associated files are retained for the lifetime of your account.
            You may delete your projects at any time. Upon account closure, data will be
            deleted within 90 days.
          </Typography>

          <Typography variant="h5" component="h2">
            Your Rights
          </Typography>
          <Typography>
            Under UK GDPR, you have the right to:
          </Typography>
          <Box component="ul" sx={{ pl: 3, mb: 2 }}>
            <li>Access your personal data</li>
            <li>Correct inaccurate data</li>
            <li>Request deletion of your data</li>
            <li>Export your data in a portable format</li>
            <li>Object to certain processing activities</li>
          </Box>

          <Typography variant="h5" component="h2">
            Contact
          </Typography>
          <Typography>
            For questions about this privacy policy or to exercise your data rights,
            please contact the Newcastle University Data Protection Officer or the
            service administrators.
          </Typography>

          <Typography variant="h5" component="h2">
            Changes to This Policy
          </Typography>
          <Typography>
            We may update this privacy policy from time to time. Significant changes
            will be communicated through the service interface.
          </Typography>
        </Box>
      </Paper>
    </Container>
  );
}
