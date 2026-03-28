/*
 * Copyright (C) 2025 Newcastle University
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
import React from "react";
import { Typography, Stack, Link } from "@mui/material";
import { CCP4i2ReportElementProps } from "./CCP4i2ReportElement";

export const CCP4i2ReportReference: React.FC<CCP4i2ReportElementProps> = ({
  item,
  iItem,
}) => {
  const title = item.getAttribute("articleTitle") || "Untitled";
  const source = item.getAttribute("source") || "Unknown source";

  let authors: string[] = [];
  const authorListRaw = item.getAttribute("authorList");

  try {
    if (authorListRaw) {
      authors = JSON.parse(authorListRaw.replace(/'/g, '"'));
    }
  } catch (e) {
    console.warn("Failed to parse author list:", authorListRaw);
  }

  // Encode the article title for the PubMed search URL
  const pubmedSearchUrl = `https://pubmed.ncbi.nlm.nih.gov/?term=${encodeURIComponent(
    title
  )}`;

  return (
    <Stack spacing={0.5} sx={{ mb: 2 }}>
      <Typography variant="body2" fontWeight="bold">
        {iItem} {authors.length ? authors.join(", ") : "Unknown authors"}
      </Typography>
      <Typography variant="body2" fontStyle="italic">
        <Link
          href={pubmedSearchUrl}
          target="_blank"
          rel="noopener noreferrer"
          underline="hover"
        >
          {title}
        </Link>
      </Typography>
      <Typography variant="body2">{source}</Typography>
    </Stack>
  );
};

export default CCP4i2ReportReference;
