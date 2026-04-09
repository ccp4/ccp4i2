'use client';

import { Container, Typography, Box, Paper, Link as MuiLink } from '@mui/material';
import Link from 'next/link';
import ArrowBackIcon from '@mui/icons-material/ArrowBack';

export default function TermsOfUsePage() {
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
          Terms of Use
        </Typography>

        <Typography variant="body2" color="text.secondary" sx={{ mb: 4 }}>
          Last updated: January 2025
        </Typography>

        <Box sx={{ '& h2': { mt: 4, mb: 2 }, '& p': { mb: 2 } }}>
          <Typography variant="h5" component="h2">
            Acceptance of Terms
          </Typography>
          <Typography>
            By accessing and using CCP4i2, you agree to be bound by these Terms of Use.
            If you do not agree to these terms, please do not use the service.
          </Typography>

          <Typography variant="h5" component="h2">
            Service Description
          </Typography>
          <Typography>
            CCP4i2 is a crystallographic computing platform that provides tools for
            macromolecular structure determination, refinement, and analysis. The service
            includes compound registration, assay management, and related scientific
            data management capabilities.
          </Typography>

          <Typography variant="h5" component="h2">
            Eligibility
          </Typography>
          <Typography>
            This service is available to authorised users with valid institutional
            credentials. Access is provided through Microsoft Entra ID authentication
            and is subject to your organisation&apos;s policies.
          </Typography>

          <Typography variant="h5" component="h2">
            Acceptable Use
          </Typography>
          <Typography>
            You agree to use this service only for legitimate scientific research purposes.
            You must not:
          </Typography>
          <Box component="ul" sx={{ pl: 3, mb: 2 }}>
            <li>Attempt to gain unauthorised access to the system or other users&apos; data</li>
            <li>Use the service to store or process illegal content</li>
            <li>Interfere with the operation of the service</li>
            <li>Share your credentials or allow others to access your account</li>
            <li>Use the service for commercial purposes without authorisation</li>
          </Box>

          <Typography variant="h5" component="h2">
            Intellectual Property
          </Typography>
          <Typography>
            <strong>Your Data:</strong> You retain ownership of all scientific data,
            structures, and compounds you upload to the service.
          </Typography>
          <Typography>
            <strong>Software:</strong> CCP4i2 incorporates the CCP4 Software Suite and
            other scientific software packages. Use of these components is subject to
            their respective licences. CCP4 is provided under the LGPL licence for
            academic use.
          </Typography>

          <Typography variant="h5" component="h2">
            Data Responsibility
          </Typography>
          <Typography>
            You are responsible for maintaining appropriate backups of your data.
            While we implement measures to protect data integrity, you should not
            rely on this service as your sole data repository for critical research data.
          </Typography>

          <Typography variant="h5" component="h2">
            Service Availability
          </Typography>
          <Typography>
            We aim to provide reliable service but do not guarantee uninterrupted
            availability. The service may be suspended for maintenance, upgrades,
            or circumstances beyond our control. We will endeavour to provide advance
            notice of planned downtime where possible.
          </Typography>

          <Typography variant="h5" component="h2">
            Limitation of Liability
          </Typography>
          <Typography>
            The service is provided &quot;as is&quot; without warranty of any kind. To the
            fullest extent permitted by law, Newcastle University shall not be liable
            for any indirect, incidental, or consequential damages arising from your
            use of the service, including loss of data or research outcomes.
          </Typography>

          <Typography variant="h5" component="h2">
            Termination
          </Typography>
          <Typography>
            We reserve the right to suspend or terminate access to the service for
            violation of these terms or for operational reasons. You may discontinue
            use at any time by contacting the service administrators.
          </Typography>

          <Typography variant="h5" component="h2">
            Changes to Terms
          </Typography>
          <Typography>
            We may modify these terms at any time. Continued use of the service after
            changes constitutes acceptance of the new terms. Significant changes will
            be communicated through the service interface.
          </Typography>

          <Typography variant="h5" component="h2">
            Governing Law
          </Typography>
          <Typography>
            These terms are governed by the laws of England and Wales. Any disputes
            shall be subject to the exclusive jurisdiction of the courts of England
            and Wales.
          </Typography>

          <Typography variant="h5" component="h2">
            Contact
          </Typography>
          <Typography>
            For questions about these terms, please contact the service administrators
            at Newcastle University.
          </Typography>
        </Box>
      </Paper>
    </Container>
  );
}
