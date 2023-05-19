import uproot as up

ufile = up.open('../output/test_At100k_160k_boxes_10.5um/out_AtDNA_100k_spacing_2um_5115.root')
eventEdep = ufile['ntuple/EventEdep;1']
edep_J = eventEdep['Edep_J'].array()
print(len(edep_J))
direct = ufile['ntuple/Direct;1']
x = direct['x'].array()
print(len(x))
print("done.")