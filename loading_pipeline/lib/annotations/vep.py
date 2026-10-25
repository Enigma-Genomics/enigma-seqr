from collections.abc import Callable

import hail as hl

from loading_pipeline.lib.annotations.enums import (
    BIOTYPES,
    FIVEUTR_CONSEQUENCES,
    LOF_FILTERS,
    TRANSCRIPT_CONSEQUENCE_TERMS,
    validated_enum_member,
)

EXTENDED_INTRONIC_SPLICE_REGION_VARIANT = 'extended_intronic_splice_region_variant'
NAGNAG_SITE = 'NAGNAG_SITE'
OMIT_TRANSCRIPT_CONSEQUENCE_TERMS = hl.set(
    [
        'upstream_gene_variant',
        'downstream_gene_variant',
    ],
)
PROTEIN_CODING_BIOTYPE = 'protein_coding'


def _lof_filters(c: hl.StructExpression) -> hl.ArrayExpression:
    return hl.or_missing(
        (c.lof == 'LC') & hl.is_defined(c.lof_filter),
        c.lof_filter.split('&|,').map(lambda f: validated_enum_member(f, LOF_FILTERS)),
    )


def _consequence_terms(c: hl.StructExpression) -> hl.ArrayExpression:
    return c.consequence_terms.filter(
        lambda t: ~OMIT_TRANSCRIPT_CONSEQUENCE_TERMS.contains(t),
    ).map(lambda t: validated_enum_member(t, TRANSCRIPT_CONSEQUENCE_TERMS))


def _legacy_utr_effects(c):
    # Retain the information old records actually contain, including every entry.
    return hl.or_else(hl.or_missing(hl.is_defined(c.fiveutr_consequence), hl.array([
        hl.struct(
            effect=hl.missing(hl.tstr),
            consequence=validated_enum_member(c.fiveutr_consequence, FIVEUTR_CONSEQUENCES),
            annotations=hl.sorted(c.fiveutr_annotation.items(), key=lambda entry: entry[0]).map(
                lambda entry: hl.struct(index=entry[0], annotation=entry[1]),
            ),
        ),
    ])), hl.empty_array(hl.tstruct(
        effect=hl.tstr, consequence=hl.tstr,
        annotations=hl.tarray(hl.tstruct(index=hl.tstr, annotation=c.fiveutr_annotation.dtype.value_type)),
    )))


def _all_utr_effects(c):
    # Sort classes and annotation entries for serialization, never severity.
    if 'fiveutr_all_effects' not in c.dtype.fields:
        return _legacy_utr_effects(c)
    effects = hl.sorted(
        hl.or_else(c.fiveutr_all_effects.items(), hl.empty_array(
            hl.ttuple(hl.tstr, c.fiveutr_all_effects.dtype.value_type),
        )), key=lambda item: item[0],
    ).map(lambda item: hl.struct(
        effect=item[0],
        consequence=validated_enum_member(item[1].consequence, FIVEUTR_CONSEQUENCES),
        annotations=hl.sorted(item[1].annotation.items(), key=lambda entry: entry[0]).map(
            lambda entry: hl.struct(index=entry[0], annotation=entry[1]),
        ),
    ))
    return hl.if_else(hl.is_defined(c.fiveutr_all_effects), effects, _legacy_utr_effects(c))


def vep_110_transcript_consequences_select(
    gencode_ensembl_to_refseq_id_mapping: hl.tdict(hl.tstr, hl.tstr),
) -> hl.StructExpression:
    return lambda c: c.select(
        aminoAcids=c.amino_acids,
        canonical=c.canonical,
        codons=c.codons,
        geneId=c.gene_id,
        hgvsc=c.hgvsc,
        hgvsp=c.hgvsp,
        transcriptId=c.transcript_id,
        maneSelect=c.mane_select,
        manePlusClinical=c.mane_plus_clinical,
        biotype=validated_enum_member(c.biotype, BIOTYPES),
        consequenceTerms=_consequence_terms(c),
        exon=hl.bind(
            lambda split: hl.or_missing(
                hl.is_defined(split),
                hl.Struct(index=split[0], total=split[1]),
            ),
            c.exon.split('/').map(hl.parse_int32),
        ),
        intron=hl.bind(
            lambda split: hl.or_missing(
                hl.is_defined(split),
                hl.Struct(index=split[0], total=split[1]),
            ),
            c.intron.split('/').map(hl.parse_int32),
        ),
        refseqTranscriptId=gencode_ensembl_to_refseq_id_mapping.get(c.transcript_id),
        alphamissense=hl.struct(
            pathogenicity=c.am_pathogenicity,
        ),
        loftee=hl.struct(
            isLofNagnag=c.lof_flags == NAGNAG_SITE,
            lofFilters=_lof_filters(c),
        ),
        spliceregion=hl.struct(
            extended_intronic_splice_region_variant=(
                hl.is_defined(c.spliceregion)
                & c.spliceregion.contains(
                    EXTENDED_INTRONIC_SPLICE_REGION_VARIANT,
                )
            ),
        ),
        utrannotator=hl.struct(
            existingInframeOorfs=c.existing_inframe_oorfs,
            existingOutofframeOorfs=c.existing_outofframe_oorfs,
            existingUorfs=c.existing_uorfs,
            fiveutrConsequences=_all_utr_effects(c).map(lambda effect: effect.consequence),
            fiveutrEffectsJson=hl.json(_all_utr_effects(c)),
            fiveutrConsequence=validated_enum_member(
                c.fiveutr_consequence,
                FIVEUTR_CONSEQUENCES,
            ),
            # Annotation documentation here:
            # https://github.com/ImperialCardioGenetics/UTRannotator?tab=readme-ov-file#the-detailed-annotation-for-each-consequence
            # NB:
            fiveutrAnnotation=c.fiveutr_annotation['1'].annotate(
                AltStopDistanceToCDS=hl.parse_int32(
                    c.fiveutr_annotation['1'].AltStopDistanceToCDS,
                ),
                CapDistanceToStart=hl.parse_int32(
                    c.fiveutr_annotation['1'].CapDistanceToStart,
                ),
                DistanceToCDS=hl.parse_int32(
                    c.fiveutr_annotation['1'].DistanceToCDS,
                ),
                DistanceToStop=hl.parse_int32(
                    c.fiveutr_annotation['1'].DistanceToStop,
                ),
                Evidence=hl.or_missing(
                    # Just in case a weird value ("NA" or anything else) propagates
                    (
                        (c.fiveutr_annotation['1'].Evidence == 'True')
                        | (c.fiveutr_annotation['1'].Evidence == 'False')
                    ),
                    hl.bool(c.fiveutr_annotation['1'].Evidence),
                ),
                StartDistanceToCDS=hl.parse_int32(
                    c.fiveutr_annotation['1'].StartDistanceToCDS,
                ),
                newSTOPDistanceToCDS=hl.parse_int32(
                    c.fiveutr_annotation['1'].newSTOPDistanceToCDS,
                ),
                alt_type_length=hl.parse_int32(
                    c.fiveutr_annotation['1'].alt_type_length,
                ),
                ref_StartDistanceToCDS=hl.parse_int32(
                    c.fiveutr_annotation['1'].ref_StartDistanceToCDS,
                ),
                ref_type_length=hl.parse_int32(
                    c.fiveutr_annotation['1'].ref_type_length,
                ),
            ),
        ),
    )


def vep_85_transcript_consequences_select(
    c: hl.StructExpression,
) -> hl.StructExpression:
    return c.select(
        aminoAcids=c.amino_acids,
        canonical=c.canonical,
        codons=c.codons,
        geneId=c.gene_id,
        hgvsc=c.hgvsc,
        hgvsp=c.hgvsp,
        transcriptId=c.transcript_id,
        biotype=validated_enum_member(c.biotype, BIOTYPES),
        consequenceTerms=_consequence_terms(c),
        loftee=hl.struct(
            isLofNagnag=c.lof_flags == NAGNAG_SITE,
            lofFilters=_lof_filters(c),
        ),
    )


def transcript_consequences_sort(
    ht: hl.Table,
) -> Callable[[hl.StructExpression], hl.StructExpression]:
    return lambda c: hl.bind(
        lambda is_coding, is_most_severe, is_canonical: (
            hl.cond(
                is_coding,
                hl.cond(
                    is_most_severe,
                    hl.cond(is_canonical, 1, 2),
                    hl.cond(is_canonical, 3, 4),
                ),
                hl.cond(
                    is_most_severe,
                    hl.cond(is_canonical, 5, 6),
                    hl.cond(is_canonical, 7, 8),
                ),
            )
        ),
        c.biotype == PROTEIN_CODING_BIOTYPE,
        hl.set(c.consequenceTerms).contains(ht.vep.most_severe_consequence),
        hl.or_else(c.canonical, 0) == 1,
    )
