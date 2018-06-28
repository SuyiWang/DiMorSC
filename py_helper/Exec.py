from subprocess import call

def DiMorSC(id=0, persist=1, dim=3):
	call(["./DiMorSC/bin/DiMorSC", str(id)+'.bin', str(id), str(persist), str(dim)])

def Triangulate(id=0, fill=0, dim=3):
	call(["./DiMorSC/bin/Triangulation", str(id), str(fill), str(dim)])