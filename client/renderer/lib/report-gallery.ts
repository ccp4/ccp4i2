import $ from "jquery";

/**
 * What a gallery lists an item as. Not its key: "Div_3" or "PhaserGraph0"
 * means nothing to a reader (#677). A graph carries its title on its data
 * (ccp4_data), not on itself, so look there before falling back to a number.
 */
export function galleryItemLabel(childItem: any, iItem: number): string {
  return (
    $(childItem).attr("title") ||
    $(childItem).find("ns0\\:ccp4_data").first().attr("title")?.trim() ||
    `Object ${iItem + 1}`
  );
}
