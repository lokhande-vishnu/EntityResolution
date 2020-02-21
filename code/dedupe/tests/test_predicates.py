import unittest
from dedupe import predicates
from future.builtins import str


class TestPuncStrip(unittest.TestCase):
    def test_sevenchar(self):
        s1 = predicates.StringPredicate(predicates.sameSevenCharStartPredicate,
                                        'foo')
        assert s1({'foo': u'fo,18v*1vaad80'}) == s1({'foo': u'fo18v1vaad80'})

    def test_set(self):
        s1 = predicates.SimplePredicate(predicates.wholeSetPredicate,
                                        'foo')
        colors = set(['red', 'blue', 'green'])
        assert s1({'foo': colors}) == (str(colors),)


class TestMetaphone(unittest.TestCase):
    def test_metaphone_token(self):
        block_val = predicates.metaphoneToken('9301 S. State St. ')
        assert block_val == set([u'STT', u'S', u'ST'])


class TestWholeSet(unittest.TestCase):
    def setUp(self):
        self.s1 = set(['red', 'blue', 'green'])

    def test_full_set(self):
        block_val = predicates.wholeSetPredicate(self.s1)
        self.assertEqual(block_val, (str(self.s1),))


class TestSetElement(unittest.TestCase):
    def setUp(self):
        self.s1 = set(['red', 'blue', 'green'])

    def test_long_set(self):
        block_val = predicates.commonSetElementPredicate(self.s1)
        self.assertEqual(set(block_val), set(('blue', 'green', 'red')))

    def test_empty_set(self):
        block_val = predicates.commonSetElementPredicate(set())
        self.assertEqual(block_val, tuple())

    def test_first_last(self):
        block_val = predicates.lastSetElementPredicate(self.s1)
        assert block_val == ('red',)
        block_val = predicates.firstSetElementPredicate(self.s1)
        assert block_val == ('blue',)

    def test_magnitude(self):
        block_val = predicates.magnitudeOfCardinality(self.s1)
        assert block_val == (u'0', )

        block_val = predicates.magnitudeOfCardinality(())
        assert block_val == ()


class TestLatLongGrid(unittest.TestCase):
    def setUp(self):
        self.latlong1 = (42.535, -5.012)

    def test_precise_latlong(self):
        block_val = predicates.latLongGridPredicate(self.latlong1)
        assert block_val == (u'[42.5, -5.0]',)
        block_val = predicates.latLongGridPredicate((0, 0))
        assert block_val == ()

    def test_exists(self):
        block_val = predicates.existsPredicate(self.latlong1)
        assert block_val == (u'1',)
        block_val = predicates.existsPredicate((0, 0))
        assert block_val == (u'0',)


class TestNumericPredicates(unittest.TestCase):
    def test_order_of_magnitude(self):
        assert predicates.orderOfMagnitude(10) == (u'1',)
        assert predicates.orderOfMagnitude(9) == (u'1',)
        assert predicates.orderOfMagnitude(2) == (u'0',)
        assert predicates.orderOfMagnitude(-2) == ()

    def test_round_to_1(self):
        assert predicates.roundTo1(22315) == (u'20000',)
        assert predicates.roundTo1(-22315) == (u'-20000',)


if __name__ == '__main__':
    unittest.main()
