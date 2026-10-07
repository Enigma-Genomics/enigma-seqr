import { getUtrEffects, getUtrConsequences } from './utrUtils'

const effects = [
  {
    effect: 'uAUG_gained',
    consequence: 'gain',
    annotations: [
      { index: '1', annotation: { Evidence: 'False' } },
      { index: '2', annotation: { DistanceToCDS: '30' } },
    ],
  },
  {
    effect: 'uSTOP_lost',
    consequence: 'loss',
    annotations: [
      { index: '1', annotation: { AltStop: 'True' } },
    ],
  },
]

describe('UTR effects', () => {
  it('retains every class and every evidence entry', () => {
    const utr = { fiveutrEffectsJson: JSON.stringify(effects) }
    expect(getUtrEffects(utr)).toEqual(effects)
    expect(getUtrConsequences(utr)).toEqual(['gain', 'loss'])
  })
  it('supports older loaded single-effect annotations', () => {
    expect(getUtrEffects({ fiveutrConsequence: 'legacy', fiveutrAnnotation: null })).toEqual([
      { effect: null, consequence: 'legacy', annotations: [{ index: '1', annotation: {} }] },
    ])
  })
  it('handles absent UTR data', () => {
    expect(getUtrEffects(null)).toEqual([])
    expect(getUtrEffects(undefined)).toEqual([])
  })
  it('does not silently discard malformed evidence', () => {
    expect(() => getUtrEffects({ fiveutrEffectsJson: 'invalid' })).toThrow()
  })
})
