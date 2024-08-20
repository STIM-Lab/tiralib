import struct
from tqdm import tqdm

class network:
    
    class node:
        def __init__(self):
            self.p = [0.0, 0.0, 0.0]        
            self.e = []

    class point:
        def __init__(self):
            self.x = 0.0
            self.y = 0.0
            self.z = 0.0
            self.r = 0.0
    
    class edge:
        def __init__(self):
            self.v = [-1, -1]
            self.p = []
            
    # load a NWT file from disk (file format deprecated)
    def load_nwt(self, filename):
        # open the binary file for reading
        f = open(filename, "rb")
        
        self.identifier = f.read(14)
        self.description = f.read(58)
        
        vertices = struct.unpack("I", f.read(4))[0]
        edges = struct.unpack("I", f.read(4))[0]
        
        # load all vertices
        V = []
        for vi in tqdm(range(vertices)):
            v = self.node()
            v.p[0] = struct.unpack("f", f.read(4))[0]
            v.p[1] = struct.unpack("f", f.read(4))[0]
            v.p[2] = struct.unpack("f", f.read(4))[0]
            Eo = struct.unpack("I", f.read(4))[0]
            Ei = struct.unpack("I", f.read(4))[0]
            connected_edges =  Eo + Ei
            for ei in range(connected_edges):
                v.e.append(int.from_bytes(f.read(4), "little"))
                
            V.append(v)
            
        # load all edges
        E = []
        for ei in tqdm(range(edges)):
            e = self.edge()
            e.v[0] = struct.unpack("I", f.read(4))[0]
            e.v[1] = struct.unpack("I", f.read(4))[0]
            
            Np = struct.unpack("I", f.read(4))[0]
            
            for pi in range(Np):
                p = self.point()
                p.x = struct.unpack("f", f.read(4))[0]
                p.y = struct.unpack("f", f.read(4))[0]
                p.z = struct.unpack("f", f.read(4))[0]
                p.r = struct.unpack("f", f.read(4))[0]
                
                e.p.append(p)
            E.append(e)
            
        
        # close the input file
        f.close()
        
        self.N = V
        self.E = E
        
    def load(self, filename):
        # open the binary file for reading
        f = open(filename, "rb")
        
        self.identifier = f.read(14)
        if self.identifier != bytearray("network ntsv  ", "utf-8"):
            print("ERROR: not an ntsv file")
            return
        self.description = f.read(58)
        
        vertices = struct.unpack("I", f.read(4))[0]
        edges = struct.unpack("I", f.read(4))[0]
        
        # load all nodes
        V = []
        for vi in tqdm(range(vertices)):
            v = self.node()
            v.p[0] = struct.unpack("f", f.read(4))[0]
            v.p[1] = struct.unpack("f", f.read(4))[0]
            v.p[2] = struct.unpack("f", f.read(4))[0]
            connected_edges = struct.unpack("I", f.read(4))[0]
            for ei in range(connected_edges):
                v.e.append(int.from_bytes(f.read(4), "little"))
                
            V.append(v)
            
        # load all edges
        E = []
        for ei in tqdm(range(edges)):
            e = self.edge()
            e.v[0] = struct.unpack("I", f.read(4))[0]
            e.v[1] = struct.unpack("I", f.read(4))[0]
            
            Np = struct.unpack("I", f.read(4))[0]
            
            for pi in range(Np):
                p = self.point()
                p.x = struct.unpack("f", f.read(4))[0]
                p.y = struct.unpack("f", f.read(4))[0]
                p.z = struct.unpack("f", f.read(4))[0]
                p.r = struct.unpack("f", f.read(4))[0]
                e.p.append(p)
            
            
            surface_bytes = struct.unpack("I", f.read(4))[0]
            SurfaceBytes_temp = f.read(surface_bytes)
            volume_bytes = struct.unpack("I", f.read(4))[0]
            VolumeBytes_temp = f.read(volume_bytes)
            E.append(e)
            
        
        # close the input file
        f.close()
        
        self.N = V
        self.E = E
        
    def save(self, filename):
        
        # open the binary file for reading
        f = open(filename, "wb")
        
        f.write(bytearray("network ntsv  ", "utf-8"))
        f.write(bytearray(self.description))
        
        f.write(struct.pack("I", len(self.N)))
        f.write(struct.pack("I", len(self.E)))
        
        # save all vertices
        for vi in tqdm(range(len(self.N))):
            v = self.N[vi]
            f.write(struct.pack("f", v.p[0]))
            f.write(struct.pack("f", v.p[1]))
            f.write(struct.pack("f", v.p[2]))
            
            connected_edges = int(len(v.e))
            f.write(struct.pack("I", connected_edges))
            for ei in range(connected_edges):
                f.write(struct.pack("I", int(v.e[ei])))
            
        # save all edges
        for ei in tqdm(range(len(self.E))):
            e = self.E[ei]
            f.write(struct.pack("I", e.v[0]))
            f.write(struct.pack("I", e.v[1]))
        
            
            Np = len(e.p)
            f.write(struct.pack("I", Np))
            
            for pi in range(Np):
                p = e.p[pi]
                f.write(struct.pack("f", p.x))
                f.write(struct.pack("f", p.y))
                f.write(struct.pack("f", p.z))
                f.write(struct.pack("f", p.r))
                
            surface_bytes = int(0)
            volume_bytes = int(0)
            f.write(struct.pack("I", surface_bytes))
            f.write(struct.pack("I", volume_bytes))
        f.close()
        
    def shape(self):
        return (len(self.N), len(self.E))
        
    
        

N = network()
N.load("test.ntsv")