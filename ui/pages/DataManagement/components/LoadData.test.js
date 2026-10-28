import { initialLoadValues } from './LoadData'

describe('loading sample type handoff', () => {
  it.each(['WGS', 'WES'])('prefills an allowed %s sample type', (sampleType) => {
    expect(initialLoadValues(`?sampleType=${sampleType}`)).toEqual({ sampleType })
  })

  it.each(['', '?sampleType=genome', '?sampleType=INVALID', '?sampleType=WGS&sampleType=WES'])('ignores an ambiguous or unsupported value: %s', (search) => {
    expect(initialLoadValues(search)).toEqual({})
  })

  it('does not allow a link to bypass validation or submit a request', () => {
    expect(initialLoadValues('?sampleType=WGS&validationsToSkip=validate_sample_type&filePath=gs://example/test.vcf.gz'))
      .toEqual({ sampleType: 'WGS' })
  })
})
