import { useEffect, useLayoutEffect, useRef, useState } from 'react';

interface UseAdaptiveTableHeightOptions {
  /**
   * Maximum height for the table when there's plenty of room
   * @default 500
   */
  maxHeight?: number;
  /**
   * Minimum height for the table (below this, double scrollbar is accepted)
   * @default 200
   */
  minHeight?: number;
  /**
   * Extra vertical space to reserve (for margins, padding, etc.)
   * @default 24
   */
  bottomPadding?: number;
}

interface UseAdaptiveTableHeightResult {
  /**
   * Ref to attach to the container element whose content precedes the table.
   * The hook measures from the top of viewport to the bottom of this ref's element.
   */
  headerRef: React.RefObject<HTMLDivElement | null>;
  /**
   * The calculated table height based on available viewport space
   */
  tableHeight: number;
  /**
   * Whether the table is currently at minimum height (double scrollbar territory)
   */
  isAtMinHeight: boolean;
}

/**
 * Hook that calculates adaptive table height based on available viewport space.
 *
 * Algorithm:
 * - If available height >= maxHeight: use maxHeight (plenty of room)
 * - If available height >= minHeight: shrink to fit (single scrollbar on table)
 * - If available height < minHeight: use minHeight (accept double scrollbar)
 *
 * @example
 * ```tsx
 * const { headerRef, tableHeight } = useAdaptiveTableHeight({ maxHeight: 600 });
 *
 * return (
 *   <Container>
 *     <Box ref={headerRef}>
 *       <PageHeader />
 *       <Typography variant="h4">Title</Typography>
 *     </Box>
 *     <Box sx={{ height: tableHeight }}>
 *       <DataTable fillHeight />
 *     </Box>
 *   </Container>
 * );
 * ```
 */
export function useAdaptiveTableHeight(
  options: UseAdaptiveTableHeightOptions = {}
): UseAdaptiveTableHeightResult {
  const {
    maxHeight = 500,
    minHeight = 200,
    bottomPadding = 24,
  } = options;

  const headerRef = useRef<HTMLDivElement>(null);
  const [tableHeight, setTableHeight] = useState(maxHeight);
  const [isAtMinHeight, setIsAtMinHeight] = useState(false);

  // Use layout effect to measure before paint
  useLayoutEffect(() => {
    const calculateHeight = () => {
      if (!headerRef.current) {
        setTableHeight(maxHeight);
        setIsAtMinHeight(false);
        return;
      }

      // Get the bottom position of the header content
      const headerRect = headerRef.current.getBoundingClientRect();
      const headerBottom = headerRect.bottom;

      // Calculate available space for the table
      const viewportHeight = window.innerHeight;
      const availableHeight = viewportHeight - headerBottom - bottomPadding;

      // Apply the algorithm
      let newHeight: number;
      let atMin = false;

      if (availableHeight >= maxHeight) {
        // Plenty of room - use maxHeight
        newHeight = maxHeight;
      } else if (availableHeight >= minHeight) {
        // Shrink to fit available space
        newHeight = availableHeight;
      } else {
        // Accept double scrollbar - use minHeight
        newHeight = minHeight;
        atMin = true;
      }

      setTableHeight(newHeight);
      setIsAtMinHeight(atMin);
    };

    // Initial calculation
    calculateHeight();

    // Create ResizeObserver to watch for changes
    const resizeObserver = new ResizeObserver(() => {
      calculateHeight();
    });

    // Observe the header element
    if (headerRef.current) {
      resizeObserver.observe(headerRef.current);
    }

    // Also listen for window resize
    window.addEventListener('resize', calculateHeight);

    return () => {
      resizeObserver.disconnect();
      window.removeEventListener('resize', calculateHeight);
    };
  }, [maxHeight, minHeight, bottomPadding]);

  return {
    headerRef,
    tableHeight,
    isAtMinHeight,
  };
}
