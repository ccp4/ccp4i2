import { describe, it, expect } from 'vitest';
import { protocolColour, protocolTier } from '../protocol-colour';

describe('protocolColour', () => {
  describe('unconfigured / invalid inputs', () => {
    it('returns null colour when value is null', () => {
      const r = protocolColour(null, { target_value: 10, poor_value: 1000 });
      expect(r.background).toBeNull();
      expect(r.t).toBeNull();
    });

    it('returns null colour when value is undefined', () => {
      const r = protocolColour(undefined, { target_value: 10, poor_value: 1000 });
      expect(r.background).toBeNull();
    });

    it('returns null colour when value is not finite', () => {
      const r = protocolColour(Infinity, { target_value: 10, poor_value: 1000 });
      expect(r.background).toBeNull();
    });

    it('returns null colour when protocol is null', () => {
      const r = protocolColour(100, null);
      expect(r.background).toBeNull();
    });

    it('returns null colour when target_value is null', () => {
      const r = protocolColour(100, { target_value: null, poor_value: 1000 });
      expect(r.background).toBeNull();
    });

    it('returns null colour when poor_value is null', () => {
      const r = protocolColour(100, { target_value: 10, poor_value: null });
      expect(r.background).toBeNull();
    });

    it('returns null colour when target equals poor (degenerate)', () => {
      const r = protocolColour(100, { target_value: 50, poor_value: 50 });
      expect(r.background).toBeNull();
    });
  });

  describe('lower-is-better (IC50-like: target 10, poor 1000)', () => {
    const proto = { target_value: 10, poor_value: 1000, threshold_scale: 'log' as const };

    it('value below target clamps to full green (t=1)', () => {
      const r = protocolColour(1, proto);
      expect(r.t).toBe(1);
      expect(r.background).toContain('hsl');
    });

    it('value at target gives t=1', () => {
      const r = protocolColour(10, proto);
      expect(r.t).toBe(1);
    });

    it('value at geometric midpoint (~100 nM) gives t=0.5', () => {
      const r = protocolColour(100, proto);
      expect(r.t).toBeCloseTo(0.5, 5);
    });

    it('value at poor gives t=0', () => {
      const r = protocolColour(1000, proto);
      expect(r.t).toBe(0);
    });

    it('value above poor clamps to full red (t=0)', () => {
      const r = protocolColour(10000, proto);
      expect(r.t).toBe(0);
    });

    it('the colour is monotonic: better values give higher t', () => {
      const r1 = protocolColour(50, proto);
      const r2 = protocolColour(200, proto);
      const r3 = protocolColour(500, proto);
      expect(r1.t!).toBeGreaterThan(r2.t!);
      expect(r2.t!).toBeGreaterThan(r3.t!);
    });
  });

  describe('higher-is-better (solubility-like: target 100, poor 10)', () => {
    const proto = { target_value: 100, poor_value: 10, threshold_scale: 'log' as const };

    it('value at target gives t=1 (green)', () => {
      const r = protocolColour(100, proto);
      expect(r.t).toBe(1);
    });

    it('value at poor gives t=0 (red)', () => {
      const r = protocolColour(10, proto);
      expect(r.t).toBe(0);
    });

    it('value at geometric midpoint (~31.6) gives t≈0.5', () => {
      const r = protocolColour(Math.sqrt(100 * 10), proto);
      expect(r.t).toBeCloseTo(0.5, 5);
    });

    it('value beyond target clamps to green', () => {
      const r = protocolColour(1000, proto);
      expect(r.t).toBe(1);
    });

    it('value beyond poor clamps to red', () => {
      const r = protocolColour(1, proto);
      expect(r.t).toBe(0);
    });

    it('the colour is monotonic: higher values give higher t', () => {
      const r1 = protocolColour(80, proto);
      const r2 = protocolColour(30, proto);
      const r3 = protocolColour(12, proto);
      expect(r1.t!).toBeGreaterThan(r2.t!);
      expect(r2.t!).toBeGreaterThan(r3.t!);
    });
  });

  describe('linear scale', () => {
    it('interpolates linearly for percent-like metrics', () => {
      const proto = { target_value: 90, poor_value: 10, threshold_scale: 'linear' as const };
      const r = protocolColour(50, proto);
      expect(r.t).toBeCloseTo(0.5, 5);
    });

    it('defaults to log when scale is unspecified', () => {
      const protoLog = { target_value: 10, poor_value: 1000 }; // no scale
      const protoExplicit = { target_value: 10, poor_value: 1000, threshold_scale: 'log' as const };
      expect(protocolColour(100, protoLog).t).toBeCloseTo(
        protocolColour(100, protoExplicit).t!,
        10,
      );
    });
  });

  describe('log scale with non-positive values', () => {
    it('falls back to linear when the value is zero or negative', () => {
      // Linear fallback: target=10, poor=1000, value=0 -> (0 - 1000) / (10 - 1000) = 1.010..., clamps to 1
      const r = protocolColour(0, { target_value: 10, poor_value: 1000, threshold_scale: 'log' });
      expect(r.t).toBe(1);
    });
  });
});

describe('protocolTier', () => {
  it('returns null when t is null', () => {
    expect(protocolTier(null)).toBeNull();
  });

  it('returns excellent when t=1 (at or beyond target)', () => {
    expect(protocolTier(1)).toBe('excellent');
    expect(protocolTier(1.0)).toBe('excellent');
  });

  it('returns good in the upper band', () => {
    expect(protocolTier(0.8)).toBe('good');
    expect(protocolTier(2 / 3)).toBe('good');
  });

  it('returns mid in the middle band', () => {
    expect(protocolTier(0.5)).toBe('mid');
    expect(protocolTier(1 / 3)).toBe('mid');
  });

  it('returns poor in the lower band', () => {
    expect(protocolTier(0.2)).toBe('poor');
  });

  it('returns failing at zero', () => {
    expect(protocolTier(0)).toBe('failing');
  });
});
