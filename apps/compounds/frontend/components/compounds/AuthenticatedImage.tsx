'use client';

import { useState, useEffect } from 'react';
import { Box, Skeleton } from '@mui/material';
import { BrokenImage } from '@mui/icons-material';

// Import auth helpers - same pattern as api.ts
let getAccessToken: () => Promise<string | null>;

try {
  const authModule = require('../../utils/auth-token');
  getAccessToken = authModule.getAccessToken;
} catch {
  getAccessToken = async () => null;
}

interface AuthenticatedImageProps {
  /** URL to fetch (must be an authenticated endpoint) */
  src: string | null | undefined;
  /** Alt text for the image */
  alt: string;
  /** Width - can be number (px) or string (e.g., '100%') */
  width?: number | string;
  /** Height - can be number (px) or string */
  height?: number | string;
  /** CSS object-fit property */
  objectFit?: 'cover' | 'contain' | 'fill' | 'none' | 'scale-down';
  /** Additional sx props for the container */
  sx?: Record<string, any>;
  /** Callback when image fails to load */
  onError?: () => void;
  /** Callback when image loads successfully */
  onLoad?: () => void;
}

/**
 * Image component that fetches images through authenticated endpoints.
 *
 * Regular <img> tags can't include auth headers, so this component:
 * 1. Fetches the image via fetch() with Authorization header
 * 2. Converts the response to a blob URL
 * 3. Renders the blob URL in an <img> tag
 *
 * Use this for images served from protected Django endpoints.
 */
export function AuthenticatedImage({
  src,
  alt,
  width = '100%',
  height = 'auto',
  objectFit = 'cover',
  sx,
  onError,
  onLoad,
}: AuthenticatedImageProps) {
  const [blobUrl, setBlobUrl] = useState<string | null>(null);
  const [loading, setLoading] = useState(true);
  const [error, setError] = useState(false);

  useEffect(() => {
    if (!src) {
      setLoading(false);
      setError(true);
      return;
    }

    let isMounted = true;
    let objectUrl: string | null = null;

    const fetchImage = async () => {
      try {
        setLoading(true);
        setError(false);

        // Get auth token
        const token = await getAccessToken();
        const headers: Record<string, string> = {};
        if (token) {
          headers['Authorization'] = `Bearer ${token}`;
        }

        // Fetch the image
        const response = await fetch(src, { headers });

        if (!response.ok) {
          throw new Error(`Failed to fetch image: ${response.status}`);
        }

        // Convert to blob
        const blob = await response.blob();

        if (!isMounted) return;

        // Create object URL
        objectUrl = URL.createObjectURL(blob);
        setBlobUrl(objectUrl);
        setLoading(false);
        onLoad?.();
      } catch (err) {
        console.error('AuthenticatedImage: Failed to load image', err);
        if (isMounted) {
          setError(true);
          setLoading(false);
          onError?.();
        }
      }
    };

    fetchImage();

    return () => {
      isMounted = false;
      if (objectUrl) {
        URL.revokeObjectURL(objectUrl);
      }
    };
  }, [src, onError, onLoad]);

  // Loading state
  if (loading) {
    return (
      <Skeleton
        variant="rectangular"
        width={width}
        height={height}
        sx={sx}
      />
    );
  }

  // Error state
  if (error || !blobUrl) {
    return (
      <Box
        sx={{
          width,
          height,
          display: 'flex',
          alignItems: 'center',
          justifyContent: 'center',
          bgcolor: 'grey.200',
          color: 'grey.500',
          ...sx,
        }}
      >
        <BrokenImage />
      </Box>
    );
  }

  // Success - render image
  return (
    <Box
      component="img"
      src={blobUrl}
      alt={alt}
      sx={{
        width,
        height,
        objectFit,
        ...sx,
      }}
    />
  );
}
