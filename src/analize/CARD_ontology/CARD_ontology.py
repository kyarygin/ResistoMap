from collections import defaultdict
import re
import os

SCRIPTDIR = os.path.dirname(os.path.realpath(__file__))

relations = [
    'confers_resistance_to',
    'confers_resistance_to_drug',
    'derives_from',
    'is_a',
    'part_of',
    'regulates',
    'targeted_by',
    'targeted_by_drug'
]

class Ontology(dict):
    def __init__(self):
        dict.__init__(self)

    def initialize(self):
        chunks = [[]]
        with open(os.path.join(SCRIPTDIR, 'data', 'aro.obo')) as f:
            for line in f:
                if line.strip():
                    chunks[-1].append(line.strip())
                else:
                    chunks.append([])
        chunks = [chunk for chunk in chunks if chunk]
        term_chunks = [chunk[1:] for chunk in chunks if chunk[0] == '[Term]']
        term_chunks = [[line.split(': ', 1) for line in chunk] for chunk in term_chunks]
        term_dict = [defaultdict(list) for chunk in term_chunks]
        for dict_, chunk in zip(term_dict, term_chunks):
            for key, value in chunk:
                dict_[key].append(value)
        for record in term_dict:
            self[record['id'][0]] = record

        for aro_id in self:
            self[aro_id]['is_a'] = [re.match('ARO:\d+', x).group(0) for x in self[aro_id]['is_a']]

        for aro_id in self:
            for relation in self[aro_id]['relationship']:
                relation_name, parent_id = re.match('([a-z_]+)\s(ARO:\d+)', relation).groups()
                if relation_name not in self[aro_id]:
                    self[aro_id][relation_name] = []
                self[aro_id][relation_name].append(parent_id)


    def get_antibiotics_names(self, head_aro_id):
        all_is_a_part_of_ids = set()
        current_aro_ids = set([head_aro_id])
        while current_aro_ids:
            all_is_a_part_of_ids.update(current_aro_ids)
            new_aro_ids = set()
            for aro_id in current_aro_ids:
                if 'is_a' in self[aro_id]:
                    new_aro_ids.update(self[aro_id]['is_a'])
                if 'part_of' in self[aro_id]:
                    new_aro_ids.update(self[aro_id]['part_of'])
            current_aro_ids = new_aro_ids.copy()

        all_confer_resistance_ids = set()
        for aro_id in all_is_a_part_of_ids:
            if 'confers_resistance_to' in self[aro_id]:
                all_confer_resistance_ids.update(self[aro_id]['confers_resistance_to'])
            if 'confers_resistance_to_drug' in self[aro_id]:
                all_confer_resistance_ids.update(self[aro_id]['confers_resistance_to_drug'])

        all_is_a_part_of_conferation_ids = set()
        current_aro_ids = all_confer_resistance_ids.copy()
        while current_aro_ids:
            all_is_a_part_of_conferation_ids.update(current_aro_ids)
            new_aro_ids = set()
            for aro_id in current_aro_ids:
                if 'is_a' in self[aro_id]:
                    new_aro_ids.update(self[aro_id]['is_a'])
                if 'part_of' in self[aro_id]:
                    new_aro_ids.update(self[aro_id]['part_of'])
            current_aro_ids = new_aro_ids.copy()

        antibiotic_molecules_ids = set(x for x in all_is_a_part_of_conferation_ids if 'is_a' in self[x] and 'ARO:1000003' in self[x]['is_a'])
        return [self[x]['name'][0] for x in antibiotic_molecules_ids]

    def gene_id_to_antibiotics(self, gene_id):
        aro_id = re.search('ARO:\d+', gene_id).group(0)
        antibiotics = self.get_antibiotics_names(aro_id)
        return antibiotics
