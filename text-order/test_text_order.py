"""Portable preprocessing and fresh ElPiGraph checks (no frozen outcomes needed)."""
from pathlib import Path
import tempfile
import unittest
import numpy as np
from threadpoolctl import threadpool_limits
import fit
import runtime as R


class TextOrderTests(unittest.TestCase):
    def test_unreadable_gap_is_not_joined_or_given_false_boundaries(self):
        meta = [dict(id='p')]
        counts, census = fit.count_csv(meta, {'p': [['ab', '?', 'cd']]}, 'characters')
        self.assertEqual(dict(counts['line']['p']), {('a', 'b'): 1, ('^', 'a'): 1, ('c', 'd'): 1, ('d', '$'): 1})
        self.assertEqual(census['readable_words']['p'], 2)
        self.assertEqual(counts['word']['p'][('$', '^')], 0)

    def test_line_boundaries_and_metadata_do_not_change_features(self):
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / 'input.csv'
            path.write_text('id,text,group\nz,ab cd,late\na,ef gh,early\nz,ab cd,late\n', encoding='utf-8')
            meta, lines = fit.read_csv(path)
            first, _ = fit.count_csv(meta, lines, 'characters')
            path.write_text('id,text,group\na,ef gh,late\nz,ab cd,early\nz,ab cd,early\n', encoding='utf-8')
            other, lines = fit.read_csv(path)
            second, _ = fit.count_csv(other, lines, 'characters')
            self.assertEqual([m['id'] for m in meta], [m['id'] for m in other])
            self.assertEqual(first, second)
            self.assertEqual(first['line']['z'][('d', 'a')], 0)

    def test_pairwise_ties_support_and_reversal(self):
        t = np.array([0., .5, .5, 1.])
        masks = [np.ones(4, bool)]
        a, _ = fit.consensus([t], masks)
        b, _ = fit.consensus([1 - t], masks)
        np.testing.assert_allclose(a + b, 1)
        self.assertEqual(a[1, 2], .5)
        with self.assertRaisesRegex(ValueError, 'jointly'):
            fit.consensus([t], [np.array([1, 1, 0, 1], bool)])
        aligned, margin = fit.align(1 - t, t)
        np.testing.assert_allclose(aligned, t)
        self.assertGreater(margin, 0)

    def test_actual_elpigraph_growth_and_projection_reversal(self):
        core = R.engine().core
        t = np.linspace(-1, 1, 14)
        x = np.column_stack([t, .2 * t ** 2])
        w = np.ones(len(x)) / len(x)
        initial, _ = fit.starts(x, w, 2, 10)
        with threadpool_limits(1):
            model, attempts = fit.best_fit(core, x, w, initial, .01, 8)
        self.assertTrue(model['converged'])
        self.assertEqual(len(model['nodes']), 8)
        self.assertEqual(len(attempts), 2)
        a = core.project_curve(x, model)['positions']
        reverse = dict(model, nodes=model['nodes'][::-1])
        b = core.project_curve(x, reverse)['positions']
        np.testing.assert_allclose(a + b, 1, atol=1e-12)
        self.assertGreater(abs(np.corrcoef(a, t)[0, 1]), .99)


if __name__ == '__main__':
    unittest.main()
