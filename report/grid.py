"""
Grid layout system for CCP4i2 reports.

This module provides MUI Grid-compatible layout containers that allow
report authors to create responsive, grid-based layouts. The grid system
maps directly to MUI Grid2 components in the React frontend.

Usage:
    from report.grid import GridContainer, GridItem, GridSpan

    # Create a two-column layout
    grid = GridContainer(spacing=2)

    # Left column (8/12 width on medium screens)
    left = GridItem(span=GridSpan(xs=12, md=8))
    left.append(table)
    left.append(graph)
    grid.append(left)

    # Right column (4/12 width on medium screens)
    right = GridItem(span=GridSpan(xs=12, md=4))
    right.append(summary)
    grid.append(right)

    # Add to report
    report.append(grid)

MUI Grid2 Mapping:
    GridContainer -> <Grid2 container spacing={} direction={}>
    GridItem -> <Grid2 size={{ xs: n, sm: n, md: n, lg: n, xl: n }}>

Breakpoints:
    xs: 0-600px (mobile)
    sm: 600-900px (tablet)
    md: 900-1200px (small desktop)
    lg: 1200-1536px (desktop)
    xl: 1536px+ (large desktop)
"""

from dataclasses import dataclass, field
from typing import Optional, List, Literal, Any
from xml.etree import ElementTree as ET
from enum import Enum

# Import from CCP4ReportParser to get all the addXxx methods (addDiv, addTable, etc.)
# This ensures GridItem can be used as a drop-in container replacement
from report.CCP4ReportParser import Container


class GridDirection(Enum):
    """Grid container flex direction."""
    ROW = "row"
    ROW_REVERSE = "row-reverse"
    COLUMN = "column"
    COLUMN_REVERSE = "column-reverse"


class GridJustify(Enum):
    """Grid container justify-content values."""
    FLEX_START = "flex-start"
    CENTER = "center"
    FLEX_END = "flex-end"
    SPACE_BETWEEN = "space-between"
    SPACE_AROUND = "space-around"
    SPACE_EVENLY = "space-evenly"


class GridAlign(Enum):
    """Grid container align-items values."""
    FLEX_START = "flex-start"
    CENTER = "center"
    FLEX_END = "flex-end"
    STRETCH = "stretch"
    BASELINE = "baseline"


@dataclass
class GridSpan:
    """
    Column span configuration for a grid item at different breakpoints.

    MUI Grid uses a 12-column system. Each attribute specifies how many
    columns the item should span at that breakpoint and above.

    Attributes:
        xs: Columns at extra-small (0-600px) - default 12 (full width)
        sm: Columns at small (600-900px) - inherits from xs if None
        md: Columns at medium (900-1200px) - inherits from sm if None
        lg: Columns at large (1200-1536px) - inherits from md if None
        xl: Columns at extra-large (1536px+) - inherits from lg if None

    Special values:
        'auto': Size based on content
        'grow': Grow to fill available space
        True: Same as 'grow'

    Examples:
        # Full width on mobile, half width on medium+
        GridSpan(xs=12, md=6)

        # Full width on mobile, third width on large+
        GridSpan(xs=12, lg=4)

        # Equal width columns (auto-sizing)
        GridSpan(xs='auto')
    """
    xs: int | str = 12
    sm: Optional[int | str] = None
    md: Optional[int | str] = None
    lg: Optional[int | str] = None
    xl: Optional[int | str] = None

    def as_dict(self) -> dict:
        """Convert to dictionary for serialization."""
        result = {'xs': self.xs}
        if self.sm is not None:
            result['sm'] = self.sm
        if self.md is not None:
            result['md'] = self.md
        if self.lg is not None:
            result['lg'] = self.lg
        if self.xl is not None:
            result['xl'] = self.xl
        return result


# Preset grid spans for common layouts
class GridPresets:
    """Common grid span presets for convenience."""

    # Full width at all sizes
    FULL = GridSpan(xs=12)

    # Half width on medium+
    HALF = GridSpan(xs=12, md=6)

    # Third width on medium+
    THIRD = GridSpan(xs=12, md=4)

    # Two-thirds width on medium+
    TWO_THIRDS = GridSpan(xs=12, md=8)

    # Quarter width on large+
    QUARTER = GridSpan(xs=12, sm=6, lg=3)

    # Auto-sizing
    AUTO = GridSpan(xs='auto')

    # Grow to fill
    GROW = GridSpan(xs='grow')


class GridContainer(Container):
    """
    MUI Grid container element.

    Creates a flex container that arranges its children in a grid layout.
    Maps to <Grid2 container> in the React frontend.

    Inherits from Container to get all addXxx methods (addDiv, addTable, etc.)

    Attributes:
        spacing: Gap between grid items (0-10, in MUI spacing units)
        direction: Flex direction (row, column, etc.)
        justify: Justify-content CSS value
        align: Align-items CSS value
        wrap: Whether to wrap items (default True)

    Example:
        # Basic two-column layout
        grid = GridContainer(spacing=2)
        grid.append(GridItem(span=GridSpan(xs=12, md=6), children=[left_content]))
        grid.append(GridItem(span=GridSpan(xs=12, md=6), children=[right_content]))

        # Centered content
        grid = GridContainer(
            spacing=3,
            justify=GridJustify.CENTER,
            align=GridAlign.CENTER
        )
    """

    def __init__(
        self,
        spacing: int = 2,
        direction: GridDirection | str = GridDirection.ROW,
        justify: GridJustify | str = GridJustify.FLEX_START,
        align: GridAlign | str = GridAlign.STRETCH,
        wrap: bool = True,
        **kwargs
    ):
        """
        Initialize a grid container.

        Args:
            spacing: Gap between items (0-10)
            direction: Flex direction
            justify: Justify-content value
            align: Align-items value
            wrap: Whether to wrap items
            **kwargs: Additional arguments passed to Container
        """
        super().__init__(**kwargs)

        self.spacing = max(0, min(10, spacing))  # Clamp to valid range

        # Handle string or enum for direction
        if isinstance(direction, str):
            direction = GridDirection(direction)
        self.direction = direction

        # Handle string or enum for justify
        if isinstance(justify, str):
            justify = GridJustify(justify)
        self.justify = justify

        # Handle string or enum for align
        if isinstance(align, str):
            align = GridAlign(align)
        self.align = align

        self.wrap = wrap

    def row(
        self,
        spacing: Optional[int] = None,
        justify: Optional[GridJustify | str] = None,
        align: Optional[GridAlign | str] = None
    ) -> 'GridRow':
        """
        Add a new row to this container.

        Convenience method for creating nested grid containers.

        Args:
            spacing: Override spacing for this row
            justify: Override justify for this row
            align: Override align for this row

        Returns:
            The created GridRow
        """
        row = GridRow(
            spacing=spacing if spacing is not None else self.spacing,
            justify=justify if justify is not None else self.justify,
            align=align if align is not None else self.align,
            parent=self
        )
        self.append(row)
        return row

    def item(
        self,
        span: Optional[GridSpan] = None,
        xs: int = 12,
        sm: Optional[int] = None,
        md: Optional[int] = None,
        lg: Optional[int] = None,
        xl: Optional[int] = None,
        **kwargs
    ) -> 'GridItem':
        """
        Add a new item to this container.

        Convenience method for creating grid items. You can either pass
        a GridSpan object or individual breakpoint values.

        Args:
            span: Column span configuration (overrides xs/sm/md/lg/xl)
            xs: Columns at extra-small (default 12)
            sm: Columns at small
            md: Columns at medium
            lg: Columns at large
            xl: Columns at extra-large
            **kwargs: Additional arguments for GridItem

        Returns:
            The created GridItem
        """
        if span is None:
            span = GridSpan(xs=xs, sm=sm, md=md, lg=lg, xl=xl)
        item = GridItem(span=span, parent=self, **kwargs)
        self.append(item)
        return item

    def as_data_etree(self) -> ET.Element:
        """Generate XML element for frontend."""
        root = super().as_data_etree()

        # Set grid container attributes
        root.set('container', 'true')
        root.set('spacing', str(self.spacing))
        root.set('direction', self.direction.value)
        root.set('justifyContent', self.justify.value)
        root.set('alignItems', self.align.value)
        root.set('wrap', 'wrap' if self.wrap else 'nowrap')

        return root


class GridRow(GridContainer):
    """
    A row in a grid layout (semantic alias for GridContainer).

    This is syntactic sugar for creating a nested GridContainer
    with row direction. Useful for complex multi-row layouts.

    Example:
        grid = GridContainer(spacing=2)

        row1 = grid.row()
        row1.item(GridPresets.HALF).append(content1)
        row1.item(GridPresets.HALF).append(content2)

        row2 = grid.row()
        row2.item(GridPresets.FULL).append(content3)
    """

    def __init__(self, **kwargs):
        # Force row direction
        kwargs['direction'] = GridDirection.ROW
        super().__init__(**kwargs)


class GridItem(Container):
    """
    MUI Grid item element.

    A single item within a grid container. The span configuration
    determines how many columns this item occupies at each breakpoint.
    Maps to <Grid2 size={{...}}> in the React frontend.

    Inherits from Container to get all addXxx methods (addDiv, addTable, etc.)

    Attributes:
        span: Column span at each breakpoint

    Example:
        # Full width on mobile, half on medium+
        item = GridItem(span=GridSpan(xs=12, md=6))
        item.append(my_table)

        # Using presets
        item = GridItem(span=GridPresets.THIRD)
        item.append(my_graph)
    """

    def __init__(
        self,
        span: Optional[GridSpan] = None,
        **kwargs
    ):
        """
        Initialize a grid item.

        Args:
            span: Column span configuration (defaults to full width)
            **kwargs: Additional arguments passed to Container
        """
        super().__init__(**kwargs)
        self.span = span or GridSpan()

    def as_data_etree(self) -> ET.Element:
        """Generate XML element for frontend."""
        root = super().as_data_etree()

        # Set grid item attributes
        root.set('item', 'true')

        # Set span values
        span_dict = self.span.as_dict()
        for key, value in span_dict.items():
            root.set(key, str(value))

        return root


# Convenience functions for common layouts

def two_column_layout(
    left_content: List[Any],
    right_content: List[Any],
    left_span: GridSpan = GridPresets.TWO_THIRDS,
    right_span: GridSpan = GridPresets.THIRD,
    spacing: int = 2
) -> GridContainer:
    """
    Create a two-column layout.

    Args:
        left_content: Elements for the left column
        right_content: Elements for the right column
        left_span: Span for left column (default 2/3)
        right_span: Span for right column (default 1/3)
        spacing: Gap between columns

    Returns:
        GridContainer with two columns
    """
    grid = GridContainer(spacing=spacing)

    left = grid.item(span=left_span)
    for element in left_content:
        left.append(element)

    right = grid.item(span=right_span)
    for element in right_content:
        right.append(element)

    return grid


def equal_columns(
    *contents: List[Any],
    min_width: int = 4,
    spacing: int = 2
) -> GridContainer:
    """
    Create equal-width columns.

    The number of columns is determined by the number of content lists.
    Columns will be full-width on mobile and divide evenly on larger screens.

    Args:
        *contents: List of elements for each column
        min_width: Minimum column width in grid units (1-12)
        spacing: Gap between columns

    Returns:
        GridContainer with equal columns
    """
    n_cols = len(contents)
    if n_cols == 0:
        return GridContainer(spacing=spacing)

    # Calculate column width (12 / n_cols, but at least min_width)
    col_width = max(min_width, 12 // n_cols)

    grid = GridContainer(spacing=spacing)
    for content in contents:
        item = grid.item(span=GridSpan(xs=12, md=col_width))
        for element in content:
            item.append(element)

    return grid


def stacked_layout(*contents: Any, spacing: int = 2) -> GridContainer:
    """
    Create a vertically stacked layout.

    Each content element gets full width and stacks vertically.

    Args:
        *contents: Elements to stack
        spacing: Gap between elements

    Returns:
        GridContainer with stacked elements
    """
    grid = GridContainer(spacing=spacing, direction=GridDirection.COLUMN)
    for content in contents:
        item = grid.item(span=GridPresets.FULL)
        item.append(content)

    return grid
