"""Reproduce the explicitly attributed AlN1977 data transcription and validate geometry."""
from pathlib import Path
import hashlib,json,warnings
from collections import Counter
from pymatgen.io.cif import CifParser
from larixite.cif_cluster import CIF_Cluster
ROOT=Path(__file__).resolve().parent
CIF_TEXT="# Derived experimental crystallographic data transcription, not a publisher CIF.\n# Source: H. Schulz and K. H. Thiemann, Solid State Communications 23 (1977)\n# 815-819. DOI: 10.1016/0038-1098(77)90959-0.\n# a,c: published p.816 Table 1, explicitly in angstroms.\n# x,y,z: p.816 Table 2 and its common-coordinate footnote.\n# x=-1/3,y=-2/3 are encoded as long decimal approximations only.\n# P63mc encodes the standard wurtzite setting described by the publication;\n# the numeric space-group identifier is not printed in the article.\n# Unit site occupancy expresses the published stoichiometric ordered model;\n# it is not a separately reported/refined occupancy measurement.\n# The isotropic U column of Table 2 is multiplied by 100 in the paper;\n# below it is divided by 100 to encode U_iso in angstrom squared.\n# Measurement temperature and pressure are unreported and left unknown.\n# No idealization, lattice rescaling, bond equalization or DFT relaxation.\n# Numeric facts transcribed with original article attribution; the article\n# itself remains subject to its publisher rights and is not redistributed.\ndata_AlN_Schulz_Thiemann_1977_transcription\nloop_\n_publ_author_name\n'Schulz, Heinz'\n'Thiemann, K. H.'\n_publ_section_title\n'Crystal structure refinement of AlN and GaN'\n_journal_name_full 'Solid State Communications'\n_journal_year 1977\n_journal_volume 23\n_journal_issue 11\n_journal_page_first 815\n_journal_page_last 819\n_journal_paper_doi '10.1016/0038-1098(77)90959-0'\n_chemical_name_common 'Aluminium nitride, wurtzite'\n_chemical_formula_sum 'Al N'\n_cell_length_a 3.110(1)\n_cell_length_b 3.110(1)\n_cell_length_c 4.980(1)\n_cell_angle_alpha 90\n_cell_angle_beta 90\n_cell_angle_gamma 120\n_cell_formula_units_Z 2\n_space_group_IT_number 186\n_space_group_name_H-M_alt 'P 63 m c'\n_space_group_name_Hall 'P 6c -2c'\n_symmetry_space_group_name_H-M 'P 63 m c'\n_diffrn_ambient_temperature ?\n_diffrn_ambient_pressure ?\n# Table 2 isotropic refinement: R=0.015, R(w)=0.021.\n# Reflection subset is not specified, so no CIF all/gt R-factor tag is assigned.\nloop_\n_atom_site_label\n_atom_site_type_symbol\n_atom_site_fract_x\n_atom_site_fract_y\n_atom_site_fract_z\n_atom_site_occupancy\n_atom_site_U_iso_or_equiv\nAl1 Al -0.3333333333333333 -0.6666666666666667 0 1 0.0038(2)\nN1 N -0.3333333333333333 -0.6666666666666667 0.3821(3) 1 0.0044(2)\n"

def main():
    path=ROOT/'source.cif'
    path.write_text(CIF_TEXT)
    t=json.loads((ROOT/'transcription.json').read_text())
    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter('always')
        strict=CifParser.from_str(CIF_TEXT,occupancy_tolerance=1.0).parse_structures(primitive=False,on_error='raise')[0]
        cluster=CIF_Cluster(ciftext=CIF_TEXT,absorber='Al')
    assert strict.is_ordered and cluster.struct.is_ordered
    assert strict.composition==cluster.struct.composition and str(strict.composition)=='Al2 N2'
    assert len(strict)==len(cluster.struct)==4
    sites=[]
    for index,(site,multiplicity,wyckoff) in enumerate(cluster.unique_sites,1):
        if site.specie.symbol!='Al':continue
        neighbors=cluster.struct.get_neighbors(site,6)
        distances=sorted(float(n.nn_distance) for n in neighbors if n.specie.symbol=='N')
        first=[d for d in distances if d<2.3]
        assert len(first)==4
        expected=[1.8891,1.8891,1.8891,1.9029]
        residuals=[a-b for a,b in zip(first,expected)]
        assert max(abs(d) for d in residuals)<0.00005
        other=sorted((float(n.nn_distance),n.specie.symbol) for n in neighbors if n.specie.symbol!='N')
        sites.append(dict(site_index_one_based=index,element='Al',fractional_coordinates=[float(x) for x in site.frac_coords],multiplicity=int(multiplicity),wyckoff=wyckoff,first_shell_element='N',first_shell_cutoff_A=2.3,coordination_number=4,distances_A=dict(Counter(f'{d:.6f}' for d in first)),full_first_shell_distance_multiset_A=first,next_same_element_distance_A=min(d for d in distances if d>=2.3),nearest_other_element_distance_A=other[0][0],nearest_other_element=other[0][1],site_weight=1.0,table3_residuals_A=residuals))
    assert len(sites)==1 and sites[0]['site_index_one_based']==1
    sha=hashlib.sha256(path.read_bytes()).hexdigest()
    provenance=dict(filename='source.cif',name='AlN wurtzite, Schulz-Thiemann1977 experimental transcription',absorber='Al',scatterer='N',edge='K',source_database='publication',source_record_id=t['source_record_id'],source_id=t['source_record_id'],doi=t['doi'],citation=t['citation'],primary_url=t['publisher_url'],primary_article_read_url=t['primary_article_read_url'],original_file_unmodified=False,derivation=t['derivation'],sha256=sha,cif_sha256=sha,recommendation='experimental publication transcription candidate; native FEFF and numerical validation pending',method='Single-crystal X-ray refinement; lattice parameters determined from25 reflection angles, atomic z from intensity refinement. Isotropic R=.015 and R(w)=.021; Table3 verifies the derived1+3 Al-N shell.',inequivalent_absorber_sites=1,absorber_site_indices=[1],measurement_temperature_K=None,measurement_pressure_GPa=None,measurement_conditions_evidence=t['conditions']['evidence'],synthesis_conditions=None,occupancy=dict(raw_values=['1','1'],interpretation='Explicit unit occupancies in this derived file encode the published stoichiometric ordered model; no occupancy refinement reported.'),license=t['rights'],transcription_record='transcription.json',notes=t['limitations'],inspection=dict(formula=str(cluster.struct.composition),ordered=True,hydrogen_present=False,site_count=4,spacegroup=cluster.struct.get_space_group_info()[0],spacegroup_number=cluster.struct.get_space_group_info()[1],absorber_sites=sites,parser_warnings=sorted(set(str(w.message) for w in caught))))
    reference=dict(id='AlN-wurtzite-Schulz1977-Al-K',name='AlN wurtzite,1977 single-crystal refinement',formula='AlN',absorber='Al',scatterer='N',edge='K',source_database='publication',source_record_id=t['source_record_id'],doi=t['doi'],cif_sha256=sha,original_file_unmodified=False,derivation=t['derivation'],absorber_site_indices=[1],shell_cutoff_A=2.3,expected_cn_by_site={'1':4},selection_basis='Complete tetrahedral Al-N shell,3x1.889119A+1x1.902858A, independently matching publishedTable3; nextN3.077142A. Published experimental parameters transcribed, not an originalCOD/publisher CIF.')
    for name,value in [('provenance.json',provenance),('reference.json',reference)]:
        (ROOT/name).write_text(json.dumps(value,indent=2,ensure_ascii=False)+'\n')
    print(json.dumps({'sha256':sha,'site_indices':[1],'first_shell_A':sites[0]['full_first_shell_distance_multiset_A'],'table3_residuals_A':sites[0]['table3_residuals_A'],'validation':'pass'},indent=2))
    return reference

if __name__=='__main__':main()
