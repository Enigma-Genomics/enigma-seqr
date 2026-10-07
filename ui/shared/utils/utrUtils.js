// Keep every effect and evidence entry; array order does not imply severity.
export const getUtrEffects = (value) => {
  const utrannotator = value || {}
  if (utrannotator.fiveutrEffectsJson) {
    return JSON.parse(utrannotator.fiveutrEffectsJson)
  }
  return utrannotator.fiveutrConsequence ? [{
    effect: null,
    consequence: utrannotator.fiveutrConsequence,
    annotations: [{ index: '1', annotation: utrannotator.fiveutrAnnotation || {} }],
  }] : []
}

export const getUtrConsequences = (utrannotator = {}) => (
  [...new Set(getUtrEffects(utrannotator).map(effect => effect.consequence))]
)
