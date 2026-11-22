import React from 'react'
import ReactDOM from 'react-dom'
import { act } from 'react-dom/test-utils'
import { Form, Field } from 'react-final-form'
import LoadData from './LoadData'

test('a changed VCF updates the read-only sample type and a missing receipt clears it', async () => {
  const wizard = LoadData.WrappedComponent({ user: {}, location: {} }).props.children[1]
  const Type = wizard.props.pages[0].fields.find(f => f.name === 'sampleType').component
  const div = document.createElement('div')
  let form
  global.fetch = jest.fn(url => Promise.resolve({
    ok: !url.includes('unknown'), status: url.includes('unknown') ? 400 : 200,
    json: () => Promise.resolve(url.includes('unknown') ? { error: 'Missing receipt' } :
      { locked: true, sampleType: url.includes('exome') ? 'WES' : 'WGS' }),
  }))
  await act(async () => {
    ReactDOM.render(<Form onSubmit={() => {}} render={(props) => {
      form = props.form
      return <div><Field name="filePath" component="input" /><Type /></div>
    }} />, div)
  })
  for (const [path, expected] of [['genome', 'WGS'], ['exome', 'WES'], ['unknown', undefined]]) {
    await act(async () => { form.change('filePath', path) })
    expect(form.getState().values.sampleType).toBe(expected)
    expect(div.textContent).toContain(expected || 'Missing receipt')
    expect(div.querySelector('button')).toBeNull()
  }
  act(() => { ReactDOM.unmountComponentAtNode(div) })
})
