import { useState, useEffect, useCallback, useRef } from 'react';

export interface ScrollDirectionResult {
  /** Current scroll direction */
  direction: 'up' | 'down' | 'idle';
  /** Current scroll position */
  scrollY: number;
  /** Whether user has scrolled past the collapse threshold */
  isPastThreshold: boolean;
  /** Whether the header should be in collapsed state */
  isCollapsed: boolean;
}

interface UseScrollDirectionOptions {
  /** Scroll threshold before collapsing begins (default: 100px) */
  threshold?: number;
  /** Minimum scroll delta to detect direction change (default: 10px) */
  minDelta?: number;
  /** Target element to observe scroll on (defaults to window) */
  target?: HTMLElement | null;
}

/**
 * Hook to detect scroll direction and manage collapsed state for sticky headers.
 *
 * @example
 * ```tsx
 * const { isCollapsed, direction } = useScrollDirection({ threshold: 100 });
 *
 * return (
 *   <Box sx={{ position: 'sticky', top: 0 }}>
 *     <Collapse in={!isCollapsed}>
 *       <FullDetailContent />
 *     </Collapse>
 *     <CompactSummary />
 *   </Box>
 * );
 * ```
 */
export function useScrollDirection(
  options: UseScrollDirectionOptions = {}
): ScrollDirectionResult {
  const {
    threshold = 100,
    minDelta = 10,
    target = null,
  } = options;

  const [direction, setDirection] = useState<'up' | 'down' | 'idle'>('idle');
  const [scrollY, setScrollY] = useState(0);
  const [isCollapsed, setIsCollapsed] = useState(false);

  const lastScrollY = useRef(0);
  const ticking = useRef(false);

  const updateScroll = useCallback(() => {
    const currentScrollY = target
      ? target.scrollTop
      : window.scrollY;

    const delta = currentScrollY - lastScrollY.current;

    // Only update direction if delta exceeds minimum threshold
    if (Math.abs(delta) >= minDelta) {
      const newDirection = delta > 0 ? 'down' : 'up';
      setDirection(newDirection);

      // Collapse logic:
      // - Always expanded when near top (within half threshold)
      // - Collapse when scrolling down past threshold
      // - Expand when scrolling up
      if (currentScrollY < threshold / 2) {
        // Near top - always expanded
        setIsCollapsed(false);
      } else if (newDirection === 'down' && currentScrollY > threshold) {
        setIsCollapsed(true);
      } else if (newDirection === 'up') {
        setIsCollapsed(false);
      }
    }

    setScrollY(currentScrollY);
    lastScrollY.current = currentScrollY;
    ticking.current = false;
  }, [target, threshold, minDelta]);

  const onScroll = useCallback(() => {
    if (!ticking.current) {
      requestAnimationFrame(updateScroll);
      ticking.current = true;
    }
  }, [updateScroll]);

  useEffect(() => {
    const scrollTarget = target || window;
    scrollTarget.addEventListener('scroll', onScroll, { passive: true });

    // Initial calculation
    updateScroll();

    return () => {
      scrollTarget.removeEventListener('scroll', onScroll);
    };
  }, [target, onScroll, updateScroll]);

  return {
    direction,
    scrollY,
    isPastThreshold: scrollY > threshold,
    isCollapsed,
  };
}
