import FormWizard from './FormWizard'
import { HttpRequestHelper } from 'shared/utils/httpRequestHelper'
jest.mock('shared/utils/httpRequestHelper', () => ({ HttpRequestHelper: jest.fn() }))
jest.mock('./FormWrapper', () => () => null)

describe('validation token propagation through real FormWizard', () => {
  let wizard; let requests; let responses
  beforeEach(() => {
    requests = []; responses = []
    HttpRequestHelper.mockImplementation((url, success) => ({ post: (value) => {
      requests.push({ url, value }); if (responses.length) success(responses.shift()); return Promise.resolve()
    } }))
    wizard = new FormWizard({ formatSubmitUrl: () => '/load' })
    wizard.setState = updater => { wizard.state = { ...wizard.state, ...(typeof updater === 'function' ? updater(wizard.state) : updater) } }
  })
  test('validation token and samples survive intermediate navigation and final submission', async () => {
    responses.push({ callsetValidation: 'signed-v1', vcfSamples: ['SYNTHETIC'] })
    await wizard.onPageSubmit('/validate')({ filePath: '/CASE.gatk.vcf.gz' })
    wizard.navigateNext(); wizard.navigateBack(); wizard.navigateNext()
    await wizard.onFormSubmit({ projects: ['synthetic-project'], filePath: '/CASE.gatk.vcf.gz' })
    expect(requests[1]).toEqual({ url: '/load', value: { callsetValidation: 'signed-v1', vcfSamples: ['SYNTHETIC'], projects: ['synthetic-project'], filePath: '/CASE.gatk.vcf.gz' } })
  })
  test('revalidation replaces prior token rather than retaining stale token', async () => {
    responses.push({ callsetValidation: 'signed-v1', vcfSamples: ['ONE'] })
    await wizard.onPageSubmit('/validate')({ filePath: '/ONE.gatk.vcf.gz' })
    wizard.navigateNext(); wizard.navigateBack()
    responses.push({ callsetValidation: 'signed-v2', vcfSamples: ['TWO'] })
    await wizard.onPageSubmit('/validate')({ filePath: '/TWO.gatk.vcf.gz' })
    await wizard.onFormSubmit({ filePath: '/TWO.gatk.vcf.gz' })
    expect(requests[2].value.callsetValidation).toBe('signed-v2')
    expect(requests[2].value.vcfSamples).toEqual(['TWO'])
  })
})
