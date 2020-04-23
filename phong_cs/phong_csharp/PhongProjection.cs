using System;
using System.IO;
using System.Collections;
using System.Collections.Generic;
using System.Runtime.InteropServices;
using g3;
using Stopwatch = System.Diagnostics.Stopwatch;


namespace Phong
{
    public enum PhongProjectionResult
    {
        Success,
        OutsideTetVolume,
        PhongProjectionFail,
        UnknownError
    };

    struct Tet
    {
        public uint[] idx;
    }

    struct Tri
    {
        public uint[] idx;
    }

    struct VectorX
    {
        public float[] data;

        public VectorX(int size)
        {
            data = new float[size];
        }

        public static VectorX operator +(VectorX v1, VectorX v2)
        {
            VectorX v = new VectorX(v1.data.Length);

            for (int i = 0; i < v1.data.Length; ++i)
                v.data[i] = v1.data[i] + v2.data[i];
            return v;
        }

        public static VectorX operator *(VectorX v1, float scale)
        {
            VectorX v = new VectorX(v1.data.Length);

            for (int i = 0; i < v1.data.Length; ++i)
                v.data[i] = v1.data[i] * scale;
            return v;
        }

        public double[] ToDouble()
        {
            var res = new double[data.Length];
            for (int i = 0; i < data.Length; ++i)
                res[i] = (double)data[i];
            return res;
        }

        public override string ToString()
        {
            return "(" + String.Join(" ", Array.ConvertAll<float, String>(data, Convert.ToString)) + ")";
        }
    }

    public class PhongProjection
    {
        const float negEps = -0.00001f;

        private IntPtr phong = IntPtr.Zero;

        [DllImport("Phong")]
        private static extern IntPtr createPhongObject(double[] V, int nV, int dim, uint[] F, int nF);

        [DllImport("Phong")]
        private static extern bool project(IntPtr phong, double[] p, int fid_start, float[] w);

        [DllImport("Phong")]
        private static extern bool deletePhongObject(IntPtr phong);

        [DllImport("Phong")]
        private static extern double getInitTime(IntPtr phong);

        double[] vertices;
        uint[] triangles;

        int triID = 0;

        // tet mesh
        VectorX[] TV;
        VectorX[] TV3D;
        Tet[] TT;
        Matrix4x4[] TMat;
        List<int>[] vertTetList;

        // surface mesh
        VectorX[] FV;
        VectorX[] FV3D;
        Tri[] FF;
        const int dim = 8;


        private int oldFid = 0;


        public void Init()
        {
            Stopwatch totalWatch = new Stopwatch();
            totalWatch.Start();

            //string tetMeshFile = Path.Combine(Application.streamingAssetsPath, ModelName + "_tet.txt");
            string triMeshFile = "C:\\Users\\arorar\\Documents\\VR_drawing_on_meshes\\Assets\\StreamingAssets\\head_tri.txt";
            //string outMeshFile = Path.Combine(Application.streamingAssetsPath, ModelName + "_out.obj");

            //string inMeshFile =
            //    LoadInsideOffsetSurface ?
            //        Path.Combine(Application.streamingAssetsPath, ModelName + "_in.obj") :
            //        "";


            if (phong != IntPtr.Zero)
            {
                deletePhongObject(phong);
            }

            Stopwatch triTetMeshWatch = new Stopwatch();
            triTetMeshWatch.Start();
            Debug.Log("Loading tet mesh...");
            //bool tetRes = LoadTetMesh(tetMeshFile);
            //Debug.Assert(tetRes);
            Debug.Log("Loading tri mesh...");
            bool triRes = LoadTriMesh(triMeshFile);
            Debug.Assert(triRes);
            triTetMeshWatch.Stop();
            
            vertices = new double[FV.Length * dim];
            triangles = new uint[FF.Length * 3];

            for (int i = 0; i < FV.Length; ++i)
                for (int j = 0; j < dim; ++j)
                    vertices[j * FV.Length + i] = (double)FV[i].data[j];

            for (int i = 0; i < FF.Length; ++i)
                for (int j = 0; j < 3; ++j)
                    triangles[j * FF.Length + i] = FF[i].idx[j];

            Debug.Log("Creating phong object...");
            Stopwatch phongCreateWatch = new Stopwatch();
            phongCreateWatch.Start();
            phong = createPhongObject(
                vertices, (int)FV.Length, dim,
                triangles, (int)FF.Length);
            phongCreateWatch.Stop();
            
            totalWatch.Stop();

            Debug.Log("Total initialization time: " + totalWatch.Elapsed.TotalSeconds);
            Debug.Log("Tri and tet mesh load time: " + triTetMeshWatch.Elapsed.TotalSeconds);
            Debug.Log("Phong Object creation time: " + phongCreateWatch.Elapsed.TotalSeconds);
            Debug.Log("Phong Init (C++) time: " + getInitTime(phong));
        }

        

        
        bool LoadTriMesh(string filename)
        {
            var data = File.ReadAllLines(filename);

            if (data.Length < 2)
                return false;

            int n = int.Parse(data[0]);
            FV = new VectorX[n];

            FV3D = new VectorX[n];
            int m = int.Parse(data[1]);
            FF = new Tri[m];

            int count = 0;

            for (int i = 2; i < 2 * n + 2; ++i)
            {
                var d = Array.ConvertAll(data[i].Split(' '), Single.Parse);
                if (i % 2 == 0)
                {
                    FV3D[count] = new VectorX[3];
                    FV3D[count].data = d;
                }
                else
                {
                    FV[count] = new VectorX(dim);
                    FV[count].data = d;
                    count++;
                }
            }

            if (count != FV.Length)
                return false;

            count = 0;

            for (int i = 2 * n + 2; i < 2 * n + m + 2; ++i)
            {
                var d = Array.ConvertAll(data[i].Split(' '), uint.Parse);
                FF[count] = new Tri();
                FF[count].idx = new uint[3];
                FF[count].idx = d;
                count = count + 1;
            }

            if (count != FF.Length)
                return false;

            return true;
        }

        public static void Main(string[] args)
        {
            Init();
        }
    }
}

