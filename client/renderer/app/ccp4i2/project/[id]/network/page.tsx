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
import { useParams } from "next/navigation";
import { ProjectNetwork } from "../../../../../components/project-network";
import { useProject } from "../../../../../utils";
import ToolBar from "../../../../../components/tool-bar";
import Container from "@mui/material/Container/Container";

export default function Page() {
  const params = useParams();
  const { id } = (params ?? {}) as { id: string };
  const { project } = useProject(parseInt(id));
  return (
    <div>
      <ToolBar />
      <Container>
        <h1>Project {project?.name} Network</h1>
        {id && <ProjectNetwork projectId={parseInt(id)} />}
      </Container>
    </div>
  );
}
