mb = """untitled
  ChemPy            3D                             0

 23 27  0  0  1  0  0  0  0  0999 V2000
   -2.0200   -5.0000   19.7780 C   0  0  0  0  0  0  0  0  0  0  0  0
   -1.2290   -4.0620   20.4210 C   0  0  0  0  0  0  0  0  0  0  0  0
    8.3460   -5.6600   16.4250 C   0  0  0  0  0  0  0  0  0  0  0  0
    2.8930   -2.9740   19.7160 C   0  0  0  0  0  0  0  0  0  0  0  0
   -1.4790   -5.8090   18.7910 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.1060   -3.9280   20.0810 C   0  0  0  0  0  0  0  0  0  0  0  0
    6.3280   -3.0300   18.8900 C   0  0  0  0  0  0  0  0  0  0  0  0
    6.6370   -6.8520   16.0870 C   0  0  0  0  0  0  0  0  0  0  0  0
    1.8580   -5.7990   17.6190 C   0  0  0  0  0  0  0  0  0  0  0  0
    4.9630   -5.0930   17.7140 C   0  0  0  0  0  0  0  0  0  0  0  0
    7.0660   -5.7770   16.7290 C   0  0  0  0  0  0  0  0  0  0  0  0
    6.3320   -4.9160   17.5360 C   0  0  0  0  0  0  0  0  0  0  0  0
    1.9030   -4.8640   18.5580 C   0  0  0  0  0  0  0  0  0  0  0  0
    2.9560   -4.0490   18.9480 C   0  0  0  0  0  0  0  0  0  0  0  0
    4.9680   -3.1560   19.1130 C   0  0  0  0  0  0  0  0  0  0  0  0
   -0.1430   -5.6740   18.4500 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.6380   -4.7300   19.0870 C   0  0  0  0  0  0  0  0  0  0  0  0
    8.6630   -6.6690   15.6230 N   0  0  0  0  0  0  0  0  0  0  0  0
    6.9630   -3.8900   18.1200 N   0  0  0  0  0  0  0  0  0  0  0  0
    4.1070   -2.4430   19.8270 N   0  0  0  0  0  0  0  0  0  0  0  0
    7.6690   -7.3660   15.4250 N   0  0  0  0  0  0  0  0  0  0  0  0
    0.6230   -6.2970   17.5610 N   0  0  0  0  0  0  0  0  0  0  0  0
    4.2710   -4.1860   18.5160 N   0  0  0  0  0  0  0  0  0  0  0  0
  1  5  4  0  0  0  0
  2  1  4  0  0  0  0
  3 18  4  0  0  0  0
  4 20  4  0  0  0  0
  5 16  4  0  0  0  0
  6  2  4  0  0  0  0
  7 19  4  0  0  0  0
  8 21  4  0  0  0  0
  9 13  4  0  0  0  0
 10 12  4  0  0  0  0
 11  3  4  0  0  0  0
 11  8  4  0  0  0  0
 12 11  1  0  0  0  0
 13 14  1  0  0  0  0
 14  4  4  0  0  0  0
 14 23  4  0  0  0  0
 15  7  4  0  0  0  0
 15 23  4  0  0  0  0
 16 22  4  0  0  0  0
 17  6  4  0  0  0  0
 17 13  4  0  0  0  0
 17 16  4  0  0  0  0
 19 12  4  0  0  0  0
 20 15  4  0  0  0  0
 21 18  4  0  0  0  0
 22  9  4  0  0  0  0
 23 10  4  0  0  0  0
M  END
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def FragIndicesToMol(oMol,indices):
    em = Chem.EditableMol(Chem.Mol())

    newIndices={}
    for i,idx in enumerate(indices):
        em.AddAtom(oMol.GetAtomWithIdx(idx))
        newIndices[idx]=i

    for i,idx in enumerate(indices):
        at = oMol.GetAtomWithIdx(idx)
        for bond in at.GetBonds():
            if bond.GetBeginAtomIdx()==idx:
                oidx = bond.GetEndAtomIdx()
            else:
                oidx = bond.GetBeginAtomIdx()
            # make sure every bond only gets added once:
            if oidx<idx:
                continue
            em.AddBond(newIndices[idx],newIndices[oidx],bond.GetBondType())
    res = em.GetMol()
    res.ClearComputedProps()
    Chem.GetSymmSSSR(res)
    res.UpdatePropertyCache(False)
    res._idxMap=newIndices
    return res

def _recursivelyModifyNs(mol,matches,indices=None):
    if indices is None:
        indices=[]
    res=None
    while len(matches) and res is None:
        tIndices=indices[:]
        nextIdx = matches.pop(0)
        tIndices.append(nextIdx)
        nm = Chem.Mol(mol.ToBinary())
        nm.GetAtomWithIdx(nextIdx).SetNoImplicit(True)
        nm.GetAtomWithIdx(nextIdx).SetNumExplicitHs(1)
        cp = Chem.Mol(nm.ToBinary())
        try:
            Chem.SanitizeMol(cp)
        except ValueError:
            res,indices = _recursivelyModifyNs(nm,matches,indices=tIndices)
        else:
            indices=tIndices
            res=cp
    return res,indices

def AdjustAromaticNs(m,nitrogenPattern='[n&D2;r5,r6]'):
    """
       default nitrogen pattern matches Ns in 5 rings and 6 rings in order to be able
       to fix: O=c1ccncc1
    """
    Chem.GetSymmSSSR(m)
    m.UpdatePropertyCache(False)

    # break non-ring bonds linking rings:
    em = Chem.EditableMol(m)
    linkers = m.GetSubstructMatches(Chem.MolFromSmarts('[r]!@[r]'))
    for a,b in linkers:
        em.RemoveBond(a,b)
    nm = em.GetMol()

    # build molecules from the fragments:
    fragLists = Chem.GetMolFrags(nm)
    frags = [FragIndicesToMol(nm,x) for x in fragLists]

    # loop through the fragments in turn and try to aromatize them:
    ok=True
    for i,frag in enumerate(frags):
        cp = Chem.Mol(frag.ToBinary())
        try:
            Chem.SanitizeMol(cp)
        except ValueError:
            matches = [x[0] for x in frag.GetSubstructMatches(Chem.MolFromSmarts(nitrogenPattern))]
            lres,indices=_recursivelyModifyNs(frag,matches)
            if not lres:
                print 'frag %d failed (%s)'%(i,str(fragLists[i]))
                ok=False
                break
            else:
                revMap={}
                for k,v in frag._idxMap.iteritems():
                    revMap[v]=k
                for idx in indices:
                    oatom = m.GetAtomWithIdx(revMap[idx])
                    oatom.SetNoImplicit(True)
                    oatom.SetNumExplicitHs(1)
    if not ok:
        return None
    return m

if __name__=='__main__':
    ms= (Chem.MolFromMolBlock(mb,False),
         Chem.MolFromSmiles('c1cnc2cc3ccnc3cc12',False),
         Chem.MolFromSmiles('c1cc2cc3ccnc3cc2n1',False),
         Chem.MolFromSmiles('O=c1ccnc(c1)-c1cnc2cc3ccnc3cc12',False),
         Chem.MolFromSmiles('O=c1ccnc(c1)-c1cc1',False),
         )
    for m in ms:
        print '#---------------------'
        try:
            m.UpdatePropertyCache(False)
            cp = Chem.Mol(m.ToBinary())
            Chem.SanitizeMol(cp)
            m = cp
            print 'fine:',Chem.MolToSmiles(m)
        except ValueError:
            print 'adjust'
            nm=AdjustAromaticNs(m)
            if nm is not None:
                Chem.SanitizeMol(nm)
                print 'fixed:',Chem.MolToSmiles(nm)
