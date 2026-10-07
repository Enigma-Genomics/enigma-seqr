"""Run the installed selection block with synthetic effects across Perl seeds."""
import json
import os
from pathlib import Path
import subprocess
import unittest

PLUGIN = Path(__file__).parents[2] / 'var/vep/plugins/UTRAnnotator.pm'


class UtrPluginSelectionTest(unittest.TestCase):
    def script(self, second='"loss"'):
        source = PLUGIN.read_text()
        block = source[source.index('  # Retain every matching class.'):source.index('  my %utr_effect;')]
        return '''use strict; use JSON::PP;
my $output_five_prime_flag; my $output_five_prime_annotation;
my %five_prime_flag=(uAUG_gained=>"gain",uSTOP_lost=>SECOND,uAUG_lost=>undef);
my %five_prime_annotation=(uAUG_gained=>{1=>{Evidence=>"False"},2=>{DistanceToCDS=>"30"}},uSTOP_lost=>{1=>{AltStop=>"True"}});
'''.replace('SECOND', second) + block + r'''
print JSON::PP->new->canonical->encode({selected=>$output_five_prime_flag,effects=>\%all_five_prime_effects});
'''

    def test_multiple_classes_and_all_evidence_repeat_across_seeds(self):
        results = [subprocess.check_output(['perl', '-e', self.script()], env={
            **os.environ, 'PERL_HASH_SEED': str(seed), 'PERL_PERTURB_KEYS': '2',
        }) for seed in range(6)]
        self.assertEqual(len(set(results)), 1)
        result = json.loads(results[0])
        self.assertIsNone(result['selected'])
        self.assertEqual(set(result['effects']), {'uAUG_gained', 'uSTOP_lost'})
        self.assertEqual(set(result['effects']['uAUG_gained']['annotation']), {'1', '2'})
        self.assertEqual(result['effects']['uAUG_gained']['annotation']['1']['Evidence'], 'False')

    def test_single_class_keeps_legacy_representative(self):
        result = json.loads(subprocess.check_output(['perl', '-e', self.script('undef')]))
        self.assertEqual(result['selected'], 'gain')
        self.assertEqual(len(result['effects']), 1)
