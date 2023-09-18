import os
import unittest

from revonto.ontology import GODag


class TestCreateGODagNoObsolete(unittest.TestCase):
    def setUp(self):
        self.godag = GODag(
            os.path.join(os.path.dirname(os.path.abspath(__file__)), "data/go1.obo")
        )

    def test_obo_dataversion(self):
        self.assertEqual(self.godag.data_version, "test/2023-09-03")

    def test_number_of_elements(self):
        self.assertEqual(len(self.godag), 7)  # 5 non-obsolete terms + 1 alt_id

    def test_GO0000006_entry(self):
        """Check all the fields of GO:0000006 entry"""
        self.assertEqual(self.godag["GO:0000006"].term_id, "GO:0000006")
        self.assertEqual(self.godag["GO:0000006"].name, "third level 2")
        self.assertEqual(self.godag["GO:0000006"].namespace, "molecular_function")
        self.assertEqual(
            self.godag["GO:0000006"].description,
            '"child of second level 1" [TC:2.A.5.1.1]',
        )
        self.assertEqual(self.godag["GO:0000006"]._parents, {"GO:0000002"})
        self.assertEqual(
            self.godag["GO:0000006"].parents, {self.godag["GO:0000002"]}
        )  # direct parents
        self.assertEqual(
            self.godag["GO:0000006"].children, {self.godag["GO:0000015"]}
        )  # direct children
        self.assertEqual(self.godag["GO:0000006"].depth, 2)
        self.assertEqual(self.godag["GO:0000006"].height, 1)

    def test_GO0000015_partof(self):
        """Check that part of relationship works for GO:0000015 which is part of GO:0005829. It also tests has_parent function"""
        self.assertEqual(
            self.godag["GO:0000015"].parents,
            {self.godag["GO:0000006"], self.godag["GO:0005829"]},
        )

    def test_has_parent(self):
        """Check has_parent function of class GOTerm"""
        self.assertTrue(self.godag["GO:0000002"].has_parent("GO:0000001"))

    def test_has_child(self):
        """Check has_child function of class GOTerm"""
        self.assertTrue(self.godag["GO:0000002"].has_child("GO:0000006"))

    def test_get_all_parents(self):
        """Check get_all_parents function of class GOTerm"""
        self.assertEqual(
            self.godag["GO:0000006"].get_all_parents(), {"GO:0000002", "GO:0000001"}
        )

    def test_get_all_children(self):
        """Check get_all_children function of class GOTerm"""
        self.assertEqual(
            self.godag["GO:0000002"].get_all_children(), {"GO:0000006", "GO:0000015"}
        )


if __name__ == "__main__":
    unittest.main()
