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
import { JobView } from "../../../../components/job-view";
import { useParams } from "next/navigation";

export default function JobPage() {
  const params = useParams();
  const { jobid } = (params ?? {}) as { jobid: string };
  return <JobView jobid={parseInt(jobid)} />;
}
