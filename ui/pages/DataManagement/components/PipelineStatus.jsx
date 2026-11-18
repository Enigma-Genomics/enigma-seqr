import React, { useEffect, useState } from 'react'
import { Container, Header, Message, Table } from 'semantic-ui-react'

const STATES = ['Waiting', 'Queued', 'Running', 'Loaded', 'Failed', 'Needs review']
const FRAME_STYLE = { width: '100%', height: '80vh', border: 0 }
const TABLE_STYLE = { marginTop: '1em', marginBottom: '2em' }

export default function PipelineStatus() {
  const [snapshot, setSnapshot] = useState(null)
  const [error, setError] = useState(false)
  useEffect(() => {
    let closed = false
    let timer
    const controller = new AbortController()
    async function refresh() {
      try {
        const res = await fetch('/api/data_management/daily_sample_status', { signal: controller.signal })
        if (!res.ok) throw new Error()
        const data = await res.json()
        if (!closed) {
          setSnapshot(data)
          setError(false)
        }
      } catch {
        if (!closed) setError(true)
      } finally {
        if (!closed) timer = setTimeout(refresh, 60000)
      }
    }
    refresh()
    return () => {
      closed = true
      controller.abort()
      clearTimeout(timer)
    }
  }, [])
  const rows = snapshot?.rows || []
  const clock = value => new Date(value).toLocaleTimeString([], { timeZone: snapshot.timezone })
  return (
    <>
      {snapshot?.enabled !== false && (
        <Container fluid style={TABLE_STYLE}>
          <Header as="h3">Samples sent today</Header>
          {(error || snapshot?.stale) && <Message warning>Daily status is delayed.</Message>}
          {!snapshot && !error && <p>Checking send history…</p>}
          {snapshot?.enabled && (
            <>
              <p>
                {snapshot.day}
                {' '}
                ·
                {' '}
                {snapshot.timezone}
                {' '}
                · Updated
                {' '}
                {clock(snapshot.updated)}
              </p>
              <strong>
                {STATES.map(s => [s, rows.filter(row => row.status === s).length])
                  .filter(([s, n]) => n || STATES.indexOf(s) < 4).map(([s, n]) => `${s} ${n}`).join(' · ')}
              </strong>
              <Table compact>
                <Table.Header>
                  <Table.Row>
                    {['Sample', 'Project', 'Type', 'Sent', 'Status', 'Pipeline run'].map(h => (
                      <Table.HeaderCell key={h}>{h}</Table.HeaderCell>
                    ))}
                  </Table.Row>
                </Table.Header>
                <Table.Body>
                  {rows.map(row => (
                    <Table.Row key={row.sample}>
                      {[row.sample, row.project, row.type, clock(row.sent), row.status, row.run || '—'].map((v, i) => (
                      // eslint-disable-next-line react/no-array-index-key
                        <Table.Cell key={i}>{v}</Table.Cell>
                      ))}
                    </Table.Row>
                  ))}
                </Table.Body>
              </Table>
              {!rows.length && <p>No samples sent today.</p>}
            </>
          )}
        </Container>
      )}
      <iframe title="Pipeline tasks" style={FRAME_STYLE} src="/luigi_ui/static/visualiser/index.html" />
    </>
  )
}
