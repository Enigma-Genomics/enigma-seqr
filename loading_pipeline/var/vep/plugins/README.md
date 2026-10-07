# Enigma UTRAnnotator correction

Vendored from the audited VEP GRCh38 container's `/plugins/UTRAnnotator.pm`.
Original SHA256: `f36bc912113acfea55846a9aeded5537c7dfdc869a49da53eccc1a20eede87bd`.
The upstream license and attribution are retained in the file.

Biological effect calculations are unchanged. The correction retains every matching
class and its complete annotation map in `5UTR_all_effects` (VEP JSON emits lowercase
`5utr_all_effects`). A single consequence/annotation is emitted only when exactly
one class matches. Array/key ordering is serialization, not clinical severity.
Full-fidelity multi-effect output requires JSON mode, as used by the seqr pipeline.
VCF legacy fields have no representative when multiple classes match.

Downstream must retain all classes and all annotation-map entries. Previously
loaded records cannot recover discarded effects without reannotation. Install this
plugin into the Enigma VEP image, then deploy a pinned image digest; application
updates alone do not replace the plugin on Dataproc.
