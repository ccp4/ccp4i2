/*
 * Copyright (C) 2026 Newcastle University
 *
 * This file is part of CCP4i2.
 *
 * CCP4i2 is free software: you can redistribute it and/or modify it
 * under the terms of the GNU Lesser General Public License version 3,
 * modified in accordance with the provisions of the license to address
 * the requirements of UK law.
 *
 * See https://www.ccp4.ac.uk/ccp4license.php for details.
 */
"use client";
import { use } from "react";
import { Container, LinearProgress } from "@mui/material";
import { JobsGrid } from "../../../../../components/jobs-grid";
import { useProject } from "../../../../../utils";

export default function JobsPage({
  params,
}: {
  params: Promise<{ id: string }>;
}) {
  const { id } = use(params);
  const { project } = useProject(parseInt(id));
  if (!project) return <LinearProgress />;
  return (
    <Container>
      <JobsGrid projectId={parseInt(id)} size={{ xs: 12, md: 4, lg: 2 }} />
    </Container>
  );
}
