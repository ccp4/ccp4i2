/*
 * Copyright (C) 2025 Newcastle University
 *
 * This file is part of CCP4i2.
 *
 * CCP4i2 is free software: you can redistribute it and/or modify it
 * under the terms of the GNU Lesser General Public License version 3,
 * modified in accordance with the provisions of the license to address
 * the requirements of UK law.
 *
 * See https://www.ccp4.ac.uk/ccp4license.php for details.
 */
import { useMemo, useRef, useState } from "react";

export function useSet<T>(iterable?: Iterable<T>) {
  const triggerRender = useState(0)[1];
  const set = useRef(new Set<T>(iterable));
  return useMemo(() => ({
    add(value) {
      if (set.current.has(value)) return;
      set.current.add(value);
      triggerRender(i => ++i);
    },
    delete(value) {
      if (!set.current.has(value)) return;
      set.current.delete(value);
      triggerRender(i => ++i);
    },
    clear() {
      if (set.current.size === 0) return;
      set.current.clear();
      triggerRender(i => ++i);
    },
    has: (value) => set.current.has(value),
    keys: () => set.current.keys(),
    values: () => set.current.values(),
    forEach: (...args) => set.current.forEach(...args),
    [Symbol.iterator]: () => set.current.values(),
    get size() { return set.current.size; },
  }) as Set<T>, []);
}
