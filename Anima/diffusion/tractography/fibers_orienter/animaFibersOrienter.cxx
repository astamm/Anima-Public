#include <tclap/CmdLine.h>

#include <animaShapesWriter.h>
#include <animaShapesReader.h>

#include <vtkSmartPointer.h>
#include <vtkPoints.h>
#include <vtkPointData.h>
#include <vtkPolyData.h>
#include <vtkGenericCell.h>
#include <vtkDoubleArray.h>

#include <itkPoolMultiThreader.h>

void OrientTracks(vtkPolyData *tracks, vtkPolyData *refTracts, vtkDoubleArray *pointsOrderArray,
                  unsigned int startIndex, unsigned int endIndex)
{
    double firstPointPosition[3];
    double secondPointPosition[3];

    vtkSmartPointer <vtkGenericCell> cell = vtkGenericCell::New();
    vtkSmartPointer <vtkGenericCell> refCell = vtkGenericCell::New();
    refTracts->GetCell(0,refCell);
    vtkPoints *refCellPts = refCell->GetPoints();
    unsigned int numRefCellPts = refCellPts->GetNumberOfPoints();

    for (unsigned int i = startIndex;i < endIndex;++i)
    {
        // Inspect i-th cell
        tracks->GetCell(i,cell);
        vtkPoints *cellPts = cell->GetPoints();
        unsigned int numCellPts = cellPts->GetNumberOfPoints();
        unsigned int minNumPts = std::min(numCellPts,numRefCellPts);
        double totalDistForward = 0;
        double totalDistBackward = 0;

        for (unsigned int j = 0;j < minNumPts;++j)
        {
            refCellPts->GetPoint(j, firstPointPosition);

            double dist = 0;
            cellPts->GetPoint(j, secondPointPosition);

            for (unsigned int l = 0;l < 3;++l)
                dist += (firstPointPosition[l] - secondPointPosition[l]) * (firstPointPosition[l] - secondPointPosition[l]);

            totalDistForward += std::sqrt(dist);
        }

        for (unsigned int j = 0;j < minNumPts;++j)
        {
            refCellPts->GetPoint(j, firstPointPosition);

            double dist = 0;
            cellPts->GetPoint(numCellPts - j, secondPointPosition);

            for (unsigned int l = 0;l < 3;++l)
                dist += (firstPointPosition[l] - secondPointPosition[l]) * (firstPointPosition[l] - secondPointPosition[l]);

            totalDistBackward += std::sqrt(dist);
        }

        if (totalDistBackward < totalDistForward)
        {
            // Need to reverse fiber
            unsigned int halfNumPts = std::floor(numCellPts / 2.0);
            for (unsigned int j = 0;j < halfNumPts;++j)
            {
                cellPts->GetPoint(j,firstPointPosition);
                cellPts->GetPoint(numCellPts - j - 1, secondPointPosition);

                cellPts->SetPoint(j,secondPointPosition);
                cellPts->SetPoint(numCellPts - j - 1, firstPointPosition);
            }
        }

        if (pointsOrderArray)
        {
            for (unsigned int j = 0;j < numCellPts;++j)
            {
                vtkIdType ptId = cell->GetPointIds()->GetId(j);
                pointsOrderArray->SetValue(ptId, j + 1);
            }
        }
    }
}

typedef struct
{
    vtkPolyData *tracks;
    vtkPolyData *refTracts;
    vtkDoubleArray *pointsOrderArray;
} ThreaderArguments;

ITK_THREAD_RETURN_FUNCTION_CALL_CONVENTION ThreadFilterer(void *arg)
{
    itk::MultiThreaderBase::WorkUnitInfo *threadArgs = (itk::MultiThreaderBase::WorkUnitInfo *)arg;
    unsigned int nbThread = threadArgs->WorkUnitID;
    unsigned int numTotalThread = threadArgs->NumberOfWorkUnits;

    ThreaderArguments *tmpArg = (ThreaderArguments *)threadArgs->UserData;
    unsigned int nbTotalCells = tmpArg->tracks->GetNumberOfCells();

    unsigned int step = nbTotalCells / numTotalThread;
    unsigned int startIndex = nbThread * step;
    unsigned int endIndex = (nbThread + 1) * step;

    if (nbThread == numTotalThread - 1)
        endIndex = nbTotalCells;

    OrientTracks(tmpArg->tracks, tmpArg->refTracts, tmpArg->pointsOrderArray, startIndex, endIndex);

    return ITK_THREAD_RETURN_DEFAULT_VALUE;
}

int main(int argc, char **argv)
{
    TCLAP::CmdLine cmd("Filters fibers from a vtp file using a label image and specifying with several -t and -f which labels should be touched or are forbidden for each fiber. INRIA / IRISA - VisAGeS/Empenn Team", ' ',ANIMA_VERSION);

    TCLAP::ValueArg<std::string> inArg("i","input","input tracks file",true,"","input tracks",cmd);
    TCLAP::ValueArg<std::string> refArg("r","ref","reference tracks name (first fiber will be used)",true,"","reference tracks",cmd);
    TCLAP::ValueArg<std::string> outArg("o","output","output tracks name",true,"","output tracks",cmd);

    TCLAP::SwitchArg pointsOrderArg("O","points-order","Output a scalar field showing points order",cmd);

    TCLAP::ValueArg<unsigned int> nbThreadsArg("T","nb-threads","Number of threads to run on (default: all available)",false,itk::MultiThreaderBase::GetGlobalDefaultNumberOfThreads(),"number of threads",cmd);
    try
    {
        cmd.parse(argc,argv);
    }
    catch (TCLAP::ArgException& e)
    {
        std::cerr << "Error: " << e.error() << "for argument " << e.argId() << std::endl;
        return EXIT_FAILURE;
    }

    anima::ShapesReader trackReader;
    trackReader.SetFileName(inArg.getValue());
    trackReader.Update();

    vtkSmartPointer <vtkPolyData> tracks = trackReader.GetOutput();

    // Get dummy cell so that it's thread safe
    vtkSmartPointer <vtkGenericCell> dummyCell = vtkGenericCell::New();
    tracks->GetCell(0,dummyCell);

    anima::ShapesReader refTrackReader;
    refTrackReader.SetFileName(refArg.getValue());
    refTrackReader.Update();

    vtkSmartPointer <vtkPolyData> refTracks = refTrackReader.GetOutput();

    vtkSmartPointer <vtkDoubleArray> pointsOrderArray = nullptr;
    if (pointsOrderArg.isSet())
    {
        vtkIdType nbTotalPts = tracks->GetNumberOfPoints();
        pointsOrderArray = vtkDoubleArray::New();
        pointsOrderArray->SetNumberOfComponents(1);
        pointsOrderArray->SetNumberOfValues(nbTotalPts);
        pointsOrderArray->SetName("Points order");
    }

    ThreaderArguments tmpStr;
    tmpStr.tracks = tracks;
    tmpStr.refTracts = refTracks;
    tmpStr.pointsOrderArray = pointsOrderArray;

    itk::PoolMultiThreader::Pointer mThreader = itk::PoolMultiThreader::New();
    mThreader->SetNumberOfWorkUnits(nbThreadsArg.getValue());
    mThreader->SetSingleMethod(ThreadFilterer,&tmpStr);
    mThreader->SingleMethodExecute();

    if (pointsOrderArg.isSet())
        tracks->GetPointData()->AddArray(pointsOrderArray);

    anima::ShapesWriter writer;
    writer.SetInputData(tracks);
    writer.SetFileName(outArg.getValue());
    std::cout << "Writing tracks: " << outArg.getValue() << std::endl;
    writer.Update();

    return EXIT_SUCCESS;
}
