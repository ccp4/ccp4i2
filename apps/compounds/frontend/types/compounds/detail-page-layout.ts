import { ReactNode } from 'react';
import { BreadcrumbItem } from '@/components/compounds/Breadcrumbs';

/**
 * Defines a single summary field shown in the collapsible header.
 * Used in the compact/collapsed view to show key information.
 */
export interface SummaryField {
  /** Label shown before the value (e.g., "MW", "Target") */
  label: string;
  /** The value to display - can be a string, number, or React node for custom rendering */
  value: ReactNode;
  /** Optional icon to show before the label */
  icon?: ReactNode;
}

/**
 * Configuration for the detail page summary shown in the sticky header.
 * This defines what appears in both collapsed and expanded states.
 */
export interface DetailSummaryConfig {
  /** Main title of the detail (e.g., "ABC-001234", "My Protocol") */
  title: string;
  /** Optional subtitle line shown below the title */
  subtitle?: string;
  /** Icon to display next to the title */
  titleIcon?: ReactNode;
  /** Summary fields shown in collapsed state (2-3 key fields) */
  fields: SummaryField[];
  /** Optional chips/tags shown after the title (e.g., status badges) */
  chips?: ReactNode;
  /** Optional action buttons (e.g., Edit, Delete) - shown in header */
  actions?: ReactNode;
}

/**
 * Props for the DetailPageLayout component.
 * This layout provides a sticky collapsible header with detail content.
 */
export interface DetailPageLayoutProps {
  /** Breadcrumb navigation items */
  breadcrumbs: BreadcrumbItem[];
  /** Summary configuration for the compact/collapsed header view */
  summary: DetailSummaryConfig;
  /** Full detail content shown when header is expanded (structure viewer, accordions, etc.) */
  detailContent: ReactNode;
  /** Optional loading state - shows skeleton in header */
  loading?: boolean;
  /** The main page content (table, etc.) - rendered below the header */
  children: ReactNode;
  /**
   * Whether to start with the header collapsed
   * @default true
   */
  defaultCollapsed?: boolean;
  /**
   * Number of summary fields to show when collapsed (on desktop)
   * @default 3
   */
  collapsedFieldCount?: number;
}
