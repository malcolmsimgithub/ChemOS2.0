import pubchempy as pcp
import re


def pubchem_lookup(smiles):
    try:
        out = pcp.get_compounds(smiles, namespace='smiles', searchtype='identity')
    except:
        return None
    if len(out) < 1:
        return None

    result = {'cid': out[0].cid, 'inchi': out[0].inchi, 'iupac_name': out[0].iupac_name}
    # Not all molecules have CAS numbers
    cas = re.search(r'\d{2,7}-\d{2}-\d', ','.join(out[0].synonyms))
    if cas:
        result['cas'] = cas.group()
        
    return result
