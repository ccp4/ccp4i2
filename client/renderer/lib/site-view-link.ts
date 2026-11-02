/**
 * The link that opens a campaign's site view, and the reading of it.
 *
 * `site` names one site, or several as `site=12,13,15` -- a site and its
 * subsites, or neighbours worth seeing at once -- viewed as one scene. The
 * order is kept: the first site gives the view its orientation and leads
 * its title (lib/campaign_scene.build_sites_scene on the server).
 */

/** The page that shows the scene of one site, or of several together. */
export function siteViewHref(campaignId: number, siteIds: number[]): string {
  return `/ccp4i2/moorhen-page/campaign/${campaignId}?summary=1&site=${siteIds.join(",")}`;
}

/** The site ids a `site` parameter names, in order, each once; [] for none. */
export function parseSiteParam(param: string | null | undefined): number[] {
  const ids: number[] = [];
  for (const token of (param ?? "").split(",")) {
    const id = Number(token.trim());
    if (token.trim() !== "" && Number.isInteger(id) && id > 0 && !ids.includes(id)) {
      ids.push(id);
    }
  }
  return ids;
}
