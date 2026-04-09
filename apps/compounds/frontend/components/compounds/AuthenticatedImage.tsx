'use client';

import { useState, useEffect, useRef } from 'react';
import { Box, Skeleton } from '@mui/material';
import { BrokenImage } from '@mui/icons-material';
import { authFetch } from '@/lib/compounds/api';

// Module-level cache for blob URLs to survive component remounts (e.g., virtual scrolling)
// Key: source URL, Value: { blobUrl, refCount }
const blobCache = new Map<string, { blobUrl: string; refCount: number }>();

// Cache cleanup: remove entries after 5 minutes of no usage
const CACHE_TTL_MS = 5 * 60 * 1000;
const cacheTimers = new Map<string, NodeJS.Timeout>();

function acquireBlobUrl(src: string, blobUrl: string): void {
  const existing = blobCache.get(src);
  if (existing) {
    existing.refCount++;
    // Clear any pending cleanup timer
    const timer = cacheTimers.get(src);
    if (timer) {
      clearTimeout(timer);
      cacheTimers.delete(src);
    }
  } else {
    blobCache.set(src, { blobUrl, refCount: 1 });
  }
}

function releaseBlobUrl(src: string): void {
  const entry = blobCache.get(src);
  if (!entry) return;

  entry.refCount--;
  if (entry.refCount <= 0) {
    // Schedule cleanup after TTL (allows remounting without refetch)
    const timer = setTimeout(() => {
      const current = blobCache.get(src);
      if (current && current.refCount <= 0) {
        URL.revokeObjectURL(current.blobUrl);
        blobCache.delete(src);
      }
      cacheTimers.delete(src);
    }, CACHE_TTL_MS);
    cacheTimers.set(src, timer);
  }
}

function getCachedBlobUrl(src: string): string | null {
  return blobCache.get(src)?.blobUrl ?? null;
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
 * 1. Fetches the image via authFetch() with Authorization header
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
  // Track which src we acquired from cache for cleanup
  const acquiredSrcRef = useRef<string | null>(null);

  useEffect(() => {
    if (!src) {
      setLoading(false);
      setError(true);
      return;
    }

    let isMounted = true;

    const fetchImage = async () => {
      try {
        // Check cache first (handles virtual scrolling remounts)
        const cached = getCachedBlobUrl(src);
        if (cached) {
          acquireBlobUrl(src, cached);
          acquiredSrcRef.current = src;
          setBlobUrl(cached);
          setLoading(false);
          onLoad?.();
          return;
        }

        setLoading(true);
        setError(false);

        // Fetch the image with auth headers
        const response = await authFetch(src);

        if (!response.ok) {
          throw new Error(`Failed to fetch image: ${response.status}`);
        }

        // Convert to blob
        const blob = await response.blob();

        if (!isMounted) return;

        // Create object URL and add to cache
        const objectUrl = URL.createObjectURL(blob);
        acquireBlobUrl(src, objectUrl);
        acquiredSrcRef.current = src;
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
      // Release our reference to the cached blob URL
      if (acquiredSrcRef.current) {
        releaseBlobUrl(acquiredSrcRef.current);
        acquiredSrcRef.current = null;
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
