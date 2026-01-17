/**
 * Campaign (Fragment Screening) API hooks.
 *
 * Provides React hooks for interacting with the ProjectGroup endpoints
 * to manage fragment screening campaigns.
 */

import useSWR, { mutate } from "swr";
import { useApi } from "../api";
import { apiPost, apiDelete, apiPatch, apiFetch } from "../api-fetch";
import {
  ProjectGroup,
  ProjectGroupDetail,
  ProjectGroupMembership,
  MemberProjectWithSummary,
  ParentFilesResponse,
  MembershipType,
} from "../types/campaigns";
import { Project } from "../types/models";

// =============================================================================
// Hook for campaign operations
// =============================================================================

export function useCampaignsApi() {
  const api = useApi();

  return {
    /**
     * List all fragment screening campaigns.
     * Filters to type=fragment_set by default.
     */
    useCampaigns(refreshInterval: number = 0) {
      return api.get<ProjectGroup[]>(
        "projectgroups?type=fragment_set",
        refreshInterval
      );
    },

    /**
     * Get a single campaign with its memberships.
     */
    useCampaign(id: number | null, refreshInterval: number = 0) {
      const endpoint = id ? `projectgroups/${id}` : null;
      return api.get<ProjectGroupDetail>(endpoint, refreshInterval);
    },

    /**
     * Get the parent project for a campaign.
     */
    useParentProject(campaignId: number | null, refreshInterval: number = 0) {
      const endpoint = campaignId
        ? `projectgroups/${campaignId}/parent_project`
        : null;
      return api.get<Project | null>(endpoint, refreshInterval);
    },

    /**
     * Get member projects with job summaries and KPIs.
     */
    useMemberProjects(campaignId: number | null, refreshInterval: number = 0) {
      const endpoint = campaignId
        ? `projectgroups/${campaignId}/member_projects`
        : null;
      return api.get<MemberProjectWithSummary[]>(endpoint, refreshInterval);
    },

    /**
     * Get parent project files (coordinates and FreeR).
     */
    useParentFiles(campaignId: number | null, refreshInterval: number = 0) {
      const endpoint = campaignId
        ? `projectgroups/${campaignId}/parent_files`
        : null;
      return api.get<ParentFilesResponse>(endpoint, refreshInterval);
    },

    // =========================================================================
    // Mutations
    // =========================================================================

    /**
     * Create a new campaign (without parent project).
     * @deprecated Use createCampaignWithParent instead
     */
    async createCampaign(name: string): Promise<ProjectGroup> {
      const result = await apiPost<ProjectGroup>("projectgroups/", {
        name,
        type: "fragment_set",
      });
      // Invalidate campaigns list
      mutate(
        (key) => typeof key === "string" && key.includes("projectgroups"),
        undefined,
        { revalidate: true }
      );
      return result;
    },

    /**
     * Create a new campaign with an auto-created parent project.
     * The parent project is created with the same name as the campaign.
     */
    async createCampaignWithParent(name: string): Promise<ProjectGroup> {
      const result = await apiPost<ProjectGroup>("projectgroups/create_with_parent/", {
        name,
        type: "fragment_set",
      });
      // Invalidate campaigns list and projects list
      mutate(
        (key) => typeof key === "string" && (key.includes("projectgroups") || key.includes("projects")),
        undefined,
        { revalidate: true }
      );
      return result;
    },

    /**
     * Update a campaign.
     */
    async updateCampaign(
      id: number,
      data: Partial<ProjectGroup>
    ): Promise<ProjectGroup> {
      const result = await apiPatch<ProjectGroup>(`projectgroups/${id}/`, data);
      // Invalidate related queries
      mutate(
        (key) => typeof key === "string" && key.includes("projectgroups"),
        undefined,
        { revalidate: true }
      );
      return result;
    },

    /**
     * Delete a campaign.
     */
    async deleteCampaign(id: number): Promise<void> {
      await apiDelete(`projectgroups/${id}/`);
      // Invalidate campaigns list
      mutate(
        (key) => typeof key === "string" && key.includes("projectgroups"),
        undefined,
        { revalidate: true }
      );
    },

    /**
     * Set the parent project for a campaign.
     * This replaces any existing parent.
     */
    async setParentProject(
      campaignId: number,
      projectId: number
    ): Promise<ProjectGroupMembership> {
      const result = await apiPost<ProjectGroupMembership>(
        `projectgroups/${campaignId}/set_parent/`,
        { project_id: projectId }
      );
      // Invalidate related queries
      mutate(
        (key) =>
          typeof key === "string" && key.includes(`projectgroups/${campaignId}`),
        undefined,
        { revalidate: true }
      );
      return result;
    },

    /**
     * Add a member project to a campaign.
     */
    async addMember(
      campaignId: number,
      projectId: number,
      type: MembershipType = "member"
    ): Promise<ProjectGroupMembership> {
      const result = await apiPost<ProjectGroupMembership>(
        `projectgroups/${campaignId}/add_member/`,
        { project_id: projectId, type }
      );
      // Invalidate related queries
      mutate(
        (key) =>
          typeof key === "string" && key.includes(`projectgroups/${campaignId}`),
        undefined,
        { revalidate: true }
      );
      return result;
    },

    /**
     * Remove a member project from a campaign.
     */
    async removeMember(campaignId: number, projectId: number): Promise<void> {
      await apiDelete(`projectgroups/${campaignId}/members/${projectId}/`);
      // Invalidate related queries
      mutate(
        (key) =>
          typeof key === "string" && key.includes(`projectgroups/${campaignId}`),
        undefined,
        { revalidate: true }
      );
    },
  };
}

// =============================================================================
// Compound SMILES lookup
// =============================================================================

/**
 * Hook to look up SMILES for a list of registry IDs.
 * Uses the compounds API if available.
 */
export function useSmilesLookup(regIds: number[]) {
  // Build query string for batch lookup
  // Use reg_number__in for django-filter lookup
  const queryString =
    regIds.length > 0 ? `reg_number__in=${regIds.join(",")}` : null;

  // This assumes compounds API is available at /api/compounds/
  // Use apiFetch for authenticated requests (required in Azure deployment)
  const { data, error, isLoading } = useSWR(
    queryString ? `/api/proxy/compounds/compounds/?${queryString}` : null,
    async (url: string) => {
      const response = await apiFetch(url);
      if (!response.ok) {
        throw new Error("Failed to fetch compound data");
      }
      const json = await response.json();
      return json.results || json;
    }
  );

  // Build lookup map
  const smilesMap: Record<number, string> = {};
  if (data && Array.isArray(data)) {
    data.forEach((compound: { reg_number?: number; smiles?: string }) => {
      if (compound.reg_number && compound.smiles) {
        smilesMap[compound.reg_number] = compound.smiles;
      }
    });
  }

  return {
    smilesMap,
    isLoading,
    error,
  };
}

// =============================================================================
// Job operations for campaigns
// =============================================================================

/**
 * Helper to create a job in a project.
 * Used for importing coordinates, FreeR, and running SubstituteLigand.
 */
export async function createJobInProject(
  projectId: number,
  taskName: string,
  title?: string
): Promise<{ new_job: { id: number; uuid: string; number: string } }> {
  return apiPost(`projects/${projectId}/create_task/`, {
    task_name: taskName,
    title,
  });
}

/**
 * Run a job via Celery.
 */
export async function runJob(
  jobId: number,
  queue: string = "default"
): Promise<void> {
  await apiPost(`jobs/${jobId}/run/`, { queue });
}
