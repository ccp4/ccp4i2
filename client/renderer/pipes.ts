export function shortDate(date: Date | string): string {
  if (typeof date === "string") date = new Date(date);
  const now = new Date();
  if (date.getFullYear() != now.getFullYear())
    return date.toLocaleString([], { year: "numeric", month: "short" });
  if (date.getMonth() != now.getMonth() || date.getDate() != now.getDate())
    return date.toLocaleString([], { month: "short", day: "numeric" });
  return date.toLocaleString([], { timeStyle: "short" });
}

export function fileSize(bytes: number): string {
  const factor = 1000;
  if (bytes < factor) return `${bytes.toPrecision(3)} Bytes`;
  let value = bytes / factor;
  for (const unit of ["KB", "MB", "GB", "TB"]) {
    if (value < factor) return `${value.toPrecision(3)} ${unit}`;
    value /= factor;
  }
  return `${value.toPrecision(3)} PB`;
}
