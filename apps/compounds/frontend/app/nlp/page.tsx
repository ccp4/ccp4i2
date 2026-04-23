'use client';

import { Box, Container } from '@mui/material';
import { PageHeader } from '@/components/compounds/PageHeader';
import { NLPPanel } from '@/components/compounds/nlp/NLPPanel';
import { routes } from '@/lib/compounds/routes';

export default function NLPPage() {
  return (
    <Box sx={{ display: 'flex', flexDirection: 'column', minHeight: '100vh' }}>
      <PageHeader
        breadcrumbs={[
          { label: 'Home', href: routes.home(), icon: 'home' },
          { label: 'Ask' },
        ]}
      />
      <Container maxWidth="lg" sx={{ py: 3, flexGrow: 1 }}>
        <NLPPanel />
      </Container>
    </Box>
  );
}
