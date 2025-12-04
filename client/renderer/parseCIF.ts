export function parseCIF(cifText: string): Record<string, any> {
  const blocks: Record<string, any> = {};
  const blockRegex = /(?:^|\n)data_(\S+)[\s\S]*?(?=(?:\n\s*data_\S+)|$)/g;
  let match: RegExpExecArray | null;

  while ((match = blockRegex.exec(cifText)) !== null) {
    const blockName = match[1];
    const blockContent = match[0]
      .replace(new RegExp(`^data_${blockName}`), "")
      .trim();
    blocks[blockName] = parseCIFBlock(blockContent);
  }
  return blocks;
}

function parseCIFBlock(blockText: string): Record<string, any> {
  const lines = blockText.split(/\r?\n/);
  const result: Record<string, any> = {};
  let i = 0;

  while (i < lines.length) {
    let line = lines[i].trim();
    if (!line || line.startsWith("#")) {
      i++;
      continue;
    }

    if (line.startsWith("loop_")) {
      i++;
      const keys: string[] = [];
      while (i < lines.length && lines[i].trim().startsWith("_")) {
        keys.push(lines[i].trim());
        i++;
      }
      const values: any[] = [];
      while (
        i < lines.length &&
        lines[i].trim() &&
        !lines[i].trim().startsWith("_") &&
        !lines[i].trim().startsWith("loop_")
      ) {
        values.push(lines[i].trim().split(/\s+/));
        i++;
      }
      result["loop_" + Object.keys(result).length] = { keys, values };
      continue;
    }

    if (line.startsWith("_")) {
      const match = line.match(/^(_\S+)\s+(.*)$/);
      if (match) {
        const key = match[1];
        let value = match[2];
        if (value.startsWith(";")) {
          value = value.slice(1) + "\n";
          i++;
          while (i < lines.length && !lines[i].startsWith(";")) {
            value += lines[i] + "\n";
            i++;
          }
          if (i < lines.length && lines[i].startsWith(";")) i++;
        }
        result[key] = value.trim();
      }
      i++;
      continue;
    }
    i++;
  }
  return result;
}
