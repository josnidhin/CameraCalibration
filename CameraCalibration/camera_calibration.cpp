#pragma once

#ifndef _DEBUG

#pragma comment(lib, "opencv_core249.lib")
#pragma comment(lib, "opencv_imgproc249.lib")
#pragma comment(lib, "opencv_highgui249.lib")
#pragma comment(lib, "opencv_ml249.lib")
#pragma comment(lib, "opencv_video249.lib")
#pragma comment(lib, "opencv_features2d249.lib")
#pragma comment(lib, "opencv_calib3d249.lib")
#pragma comment(lib, "opencv_objdetect249.lib")
#pragma comment(lib, "opencv_contrib249.lib")
#pragma comment(lib, "opencv_legacy249.lib")
#pragma comment(lib, "opencv_flann249.lib")

#else

#pragma comment(lib, "opencv_core249d.lib")
#pragma comment(lib, "opencv_imgproc249d.lib")
#pragma comment(lib, "opencv_highgui249d.lib")
#pragma comment(lib, "opencv_ml249d.lib")
#pragma comment(lib, "opencv_video249d.lib")
#pragma comment(lib, "opencv_features2d249d.lib")
#pragma comment(lib, "opencv_calib3d249d.lib")
#pragma comment(lib, "opencv_objdetect249d.lib")
#pragma comment(lib, "opencv_contrib249d.lib")
#pragma comment(lib, "opencv_legacy249d.lib")
#pragma comment(lib, "opencv_flann249d.lib")

#endif // !_DEBUG

#include <iostream>
#include <sstream>
#include <stdio.h>
#include <time.h>

#include <opencv2/core/core.hpp>
#include <opencv2/imgproc/imgproc.hpp>
#include <opencv2/calib3d/calib3d.hpp>
#include <opencv2/highgui/highgui.hpp>

#ifndef _CRT_SECURE_NO_WARNINGS
# define _CRT_SECURE_NO_WARNINGS
#endif

using namespace cv;
using namespace std;

enum { DETECTION = 0, CAPTURING = 1, CALIBRATED = 2 };

bool printIntrinsic = false,
	printDebugPoints = false;

class Settings
{
public:
	Settings() : goodInput(false) {}
	
	enum Pattern { NOT_EXISTING, CHESSBOARD, CIRCLES_GRID, ASYMMETRIC_CIRCLES_GRID };
	enum InputType { INVALID, CAMERA, VIDEO_FILE, IMAGE_LIST };

	/**
	 *
	 */
	void write(FileStorage& fs) const                        //Write serialization for this class
	{
		fs << "{" << "BoardSize_Width" << boardSize.width
			<< "BoardSize_Height" << boardSize.height
			<< "Square_Size" << squareSize
			//<< "Calibrate_Pattern" << patternToUse
			//<< "Calibrate_NrOfFrameToUse" << nrFrames
			//<< "Calibrate_FixAspectRatio" << aspectRatio
			//<< "Calibrate_AssumeZeroTangentialDistortion" << calibZeroTangentDist
			//<< "Calibrate_FixPrincipalPointAtTheCenter" << calibFixPrincipalPoint

			//<< "Write_DetectedFeaturePoints" << bwritePoints
			//<< "Write_extrinsicParameters" << bwriteExtrinsics
			<< "Write_outputFileName" << outputFileName

			//<< "Show_UndistortedImage" << showUndistorsed

			//<< "Input_FlipAroundHorizontalAxis" << flipVertical
			//<< "Input_Delay" << delay
			<< "Input" << input
			<< "}";
	}
	
	/**
	 *
	 */
	void read(const FileNode& node)                          //Read serialization for this class
	{
		node["BoardSize_Width"] >> boardSize.width;
		node["BoardSize_Height"] >> boardSize.height;
		//node["Calibrate_Pattern"] >> patternToUse;
		node["Square_Size"] >> squareSize;
		//node["Calibrate_NrOfFrameToUse"] >> nrFrames;
		//node["Calibrate_FixAspectRatio"] >> aspectRatio;
		//node["Write_DetectedFeaturePoints"] >> bwritePoints;
		//node["Write_extrinsicParameters"] >> bwriteExtrinsics;
		node["Write_outputFileName"] >> outputFileName;
		//node["Calibrate_AssumeZeroTangentialDistortion"] >> calibZeroTangentDist;
		//node["Calibrate_FixPrincipalPointAtTheCenter"] >> calibFixPrincipalPoint;
		//node["Input_FlipAroundHorizontalAxis"] >> flipVertical;
		//node["Show_UndistortedImage"] >> showUndistorsed;
		node["Input"] >> input;
		//node["Input_Delay"] >> delay;
		interprate();
	}
	
	/**
	 *
	 */
	void interprate()
	{
		goodInput = true;
		if (boardSize.width <= 0 || boardSize.height <= 0)
		{
			cerr << "Invalid Board size: " << boardSize.width << " " << boardSize.height << endl;
			goodInput = false;
		}
		if (squareSize <= 10e-6)
		{
			cerr << "Invalid square size " << squareSize << endl;
			goodInput = false;
		}

		if (input.empty())      // Check for valid input
			inputType = INVALID;
		else
		{
			if (readStringList(input, imageList))
			{
				inputType = IMAGE_LIST;
				nrFrames = (nrFrames < (int)imageList.size()) ? nrFrames : (int)imageList.size();
			}
		}

		if (inputType == INVALID)
		{
			cerr << " Inexistent input: " << input;
			goodInput = false;
		}

		flag = 0;
		flag |= CV_CALIB_FIX_PRINCIPAL_POINT;
		flag |= CV_CALIB_ZERO_TANGENT_DIST;
		//flag |= CV_CALIB_FIX_ASPECT_RATIO;

		calibrationPattern = CHESSBOARD;
		atImageList = 0;

	}

	/**
	 *
	 */
	Mat nextImage()
	{
		Mat result;
		if (atImageList < (int)imageList.size())
			result = imread(imageList[atImageList++], CV_LOAD_IMAGE_COLOR);

		return result;
	}

	static bool readStringList(const string& filename, vector<string>& l)
	{
		l.clear();
		FileStorage fs(filename, FileStorage::READ);
		if (!fs.isOpened())
			return false;
		FileNode n = fs.getFirstTopLevelNode();
		if (n.type() != FileNode::SEQ)
			return false;
		FileNodeIterator it = n.begin(), it_end = n.end();
		for (; it != it_end; ++it)
			l.push_back((string)*it);
		return true;
	}
public:
	Size boardSize;            // The size of the board -> Number of items by width and height
	Pattern calibrationPattern;// One of the Chessboard, circles, or asymmetric circle pattern
	float squareSize;          // The size of a square in your defined unit (point, millimeter,etc).
	//float aspectRatio;         // The aspect ratio
	//int delay;                 // In case of a video input
	//bool bwritePoints;         //  Write detected feature points
	//bool bwriteExtrinsics;     // Write extrinsic parameters
	//bool calibZeroTangentDist; // Assume zero tangential distortion
	//bool calibFixPrincipalPoint;// Fix the principal point at the center
	//bool flipVertical;          // Flip the captured images around the horizontal axis
	string outputFileName;      // The name of the file where to write
	//bool showUndistorsed;       // Show undistorted images after calibration
	string input;               // The input ->

	int cameraID;
	vector<string> imageList;
	int atImageList;
	//VideoCapture inputCapture;
	InputType inputType;
	bool goodInput;
	int flag;

private:
	int nrFrames;              // The number of frames to use from the input for calibration
	//string patternToUse;
};

/**
 *
 */
static void PrintHelp()
{
	cout << "Usage: CameraCalibration <config-file>" << endl
		<< " -c \t Create Config File (Follow on screen instruction)" << endl
		<< " -d \t Print debug info" << endl
		<< " -i \t Print intrinsic parameters" << endl;
		//<< "Near the sample file you'll find the configuration file, which has detailed help of "
		//"how to edit it.  It may be any OpenCV supported file format XML/YAML." << endl;
}

/**
 *
 */
static void read(const FileNode& node, Settings& x, const Settings& default_value = Settings())
{
	if (node.empty())
		x = default_value;
	else
		x.read(node);
}

/**
 *
 */
static double computeReprojectionErrors(const vector<vector<Point3f> >& objectPoints,
	const vector<vector<Point2f> >& imagePoints,
	const vector<Mat>& rvecs, const vector<Mat>& tvecs,
	const Mat& cameraMatrix, const Mat& distCoeffs,
	vector<float>& perViewErrors)
{
	vector<Point2f> imagePoints2;
	int i, totalPoints = 0;
	double totalErr = 0, err;
	perViewErrors.resize(objectPoints.size());

	for (i = 0; i < (int)objectPoints.size(); ++i)
	{
		projectPoints(Mat(objectPoints[i]), rvecs[i], tvecs[i], cameraMatrix,
			distCoeffs, imagePoints2);
		err = norm(Mat(imagePoints[i]), Mat(imagePoints2), CV_L2);

		int n = (int)objectPoints[i].size();
		perViewErrors[i] = (float)std::sqrt(err*err / n);
		totalErr += err*err;
		totalPoints += n;
	}

	return std::sqrt(totalErr / totalPoints);
}

/**
 *
 */
static void calcBoardCornerPositions(Size boardSize, float squareSize, vector<Point3f>& corners,
	Settings::Pattern patternType /*= Settings::CHESSBOARD*/)
{
	corners.clear();

	switch (patternType)
	{
	case Settings::CHESSBOARD:
	case Settings::CIRCLES_GRID:
		for (int i = 0; i < boardSize.height; ++i)
			for (int j = 0; j < boardSize.width; ++j)
				corners.push_back(Point3f(float(j*squareSize), float(i*squareSize), 0));
		break;

	case Settings::ASYMMETRIC_CIRCLES_GRID:
		for (int i = 0; i < boardSize.height; i++)
			for (int j = 0; j < boardSize.width; j++)
				corners.push_back(Point3f(float((2 * j + i % 2)*squareSize), float(i*squareSize), 0));
		break;
	default:
		break;
	}
}

/**
 *
 */
static bool runCalibration(Settings& s, Size& imageSize, Mat& cameraMatrix, Mat& distCoeffs,
	vector<vector<Point2f> > imagePoints, vector<Mat>& rvecs, vector<Mat>& tvecs,
	vector<float>& reprojErrs, double& totalAvgErr)
{

	cameraMatrix = Mat::eye(3, 3, CV_64F);
	if (s.flag & CV_CALIB_FIX_ASPECT_RATIO)
		cameraMatrix.at<double>(0, 0) = 1.0;

	distCoeffs = Mat::zeros(8, 1, CV_64F);

	vector<vector<Point3f> > objectPoints(1);
	calcBoardCornerPositions(s.boardSize, s.squareSize, objectPoints[0], s.calibrationPattern);

	objectPoints.resize(imagePoints.size(), objectPoints[0]);

	//Find intrinsic and extrinsic camera parameters
	double rms = calibrateCamera(objectPoints, imagePoints, imageSize, cameraMatrix,
		distCoeffs, rvecs, tvecs, s.flag | CV_CALIB_FIX_K4 | CV_CALIB_FIX_K5);

	cout << "Re-projection error reported by calibrateCamera: " << rms << endl;

	bool ok = checkRange(cameraMatrix) && checkRange(distCoeffs);

	totalAvgErr = computeReprojectionErrors(objectPoints, imagePoints,
		rvecs, tvecs, cameraMatrix, distCoeffs, reprojErrs);

	return ok;
}

/**
 * Print camera parameters to the output file
 */
static void saveCameraParams(Settings& s, Size& imageSize, Mat& cameraMatrix, Mat& distCoeffs,
	const vector<Mat>& rvecs, const vector<Mat>& tvecs,
	const vector<float>& reprojErrs, const vector<vector<Point2f> >& imagePoints,
	double totalAvgErr)
{
	FileStorage fs(s.outputFileName, FileStorage::WRITE);

	time_t tm;
	time(&tm);
	struct tm *t2 = localtime(&tm);
	char buf[1024];
	strftime(buf, sizeof(buf) - 1, "%c", t2);

	fs << "calibration_Time" << buf;

	if (!rvecs.empty() || !reprojErrs.empty())
		fs << "nrOfFrames" << (int)std::max(rvecs.size(), reprojErrs.size());
	fs << "image_Width" << imageSize.width;
	fs << "image_Height" << imageSize.height;
	fs << "board_Width" << s.boardSize.width;
	fs << "board_Height" << s.boardSize.height;
	fs << "square_Size" << s.squareSize;

	//if (s.flag & CV_CALIB_FIX_ASPECT_RATIO)
	//	fs << "FixAspectRatio" << s.aspectRatio;

	if (s.flag)
	{
		sprintf(buf, "flags: %s%s%s%s",
			s.flag & CV_CALIB_USE_INTRINSIC_GUESS ? " +use_intrinsic_guess" : "",
			s.flag & CV_CALIB_FIX_ASPECT_RATIO ? " +fix_aspectRatio" : "",
			s.flag & CV_CALIB_FIX_PRINCIPAL_POINT ? " +fix_principal_point" : "",
			s.flag & CV_CALIB_ZERO_TANGENT_DIST ? " +zero_tangent_dist" : "");
		cvWriteComment(*fs, buf, 0);

	}

	fs << "flagValue" << s.flag;

	fs << "Camera_Matrix" << cameraMatrix;
	fs << "Distortion_Coefficients" << distCoeffs;

	fs << "Avg_Reprojection_Error" << totalAvgErr;
	if (!reprojErrs.empty())
		fs << "Per_View_Reprojection_Errors" << Mat(reprojErrs);

	if (!rvecs.empty() && !tvecs.empty())
	{
		CV_Assert(rvecs[0].type() == tvecs[0].type());
		Mat bigmat((int)rvecs.size(), 6, rvecs[0].type());
		for (int i = 0; i < (int)rvecs.size(); i++)
		{
			Mat r = bigmat(Range(i, i + 1), Range(0, 3));
			Mat t = bigmat(Range(i, i + 1), Range(3, 6));

			CV_Assert(rvecs[i].rows == 3 && rvecs[i].cols == 1);
			CV_Assert(tvecs[i].rows == 3 && tvecs[i].cols == 1);
			//*.t() is MatExpr (not Mat) so we can use assignment operator
			r = rvecs[i].t();
			t = tvecs[i].t();
		}
		cvWriteComment(*fs, "a set of 6-tuples (rotation vector + translation vector) for each view", 0);
		fs << "Extrinsic_Parameters" << bigmat;
	}

	if (!imagePoints.empty())
	{
		Mat imagePtMat((int)imagePoints.size(), (int)imagePoints[0].size(), CV_32FC2);
		for (int i = 0; i < (int)imagePoints.size(); i++)
		{
			Mat r = imagePtMat.row(i).reshape(2, imagePtMat.cols);
			Mat imgpti(imagePoints[i]);
			imgpti.copyTo(r);
		}
		fs << "Image_points" << imagePtMat;
	}
}

/**
 *
 */
static bool runCalibrationAndSave(Settings& s, Size imageSize, Mat&  cameraMatrix, Mat& distCoeffs, vector<vector<Point2f> > imagePoints)
{
	vector<Mat> rvecs, tvecs;
	vector<float> reprojErrs;
	double totalAvgErr = 0;

	bool ok = runCalibration(s, imageSize, cameraMatrix, distCoeffs, imagePoints, rvecs, tvecs, reprojErrs, totalAvgErr);

	cout << (ok ? "Calibration succeeded" : "Calibration failed")
		<< ". avg re projection error = " << totalAvgErr << endl;

	if (ok)
	{
		if (printIntrinsic)
		{
			cout << "Camera Matrix" << endl;
			cout << format(cameraMatrix, "CSV") << endl;
			cout << "Distortion Coefficient" << endl;
			cout << format(distCoeffs, "CSV") << endl;
		}

		saveCameraParams(s, imageSize, cameraMatrix, distCoeffs, rvecs, tvecs, reprojErrs, imagePoints, totalAvgErr);
	}
	return ok;
}

/**
 *
 */
static void CreateConfigFile()
{
	string fileName, imgInputFile, paramOutputFile;
	int width, height, squareSize;

	cout << "Please enter a file name:" << endl;
	cin >> fileName;
	fileName += ".xml";

	cout << "Enter chess board width:" << endl;
	cin >> width;
	
	cout << "Enter chess board height:" << endl;
	cin >> height;

	cout << "Enter chess board square size:" << endl;
	cin >> squareSize;

	cout << "Enter image input file name:" << endl;
	cin >> imgInputFile;

	cout << "Enter parameter output file name:" << endl;
	cin >> paramOutputFile;

	FileStorage fs(fileName, FileStorage::WRITE);

	fs << "Settings" << "{";
	fs << "BoardSize_Width" << width;
	fs << "BoardSize_Height" << height;
	fs << "Square_Size" << squareSize;
	fs << "Input" << imgInputFile;
	fs << "Write_outputFileName" << paramOutputFile;
	fs << "}";

	fs.release();

}

/**
 *
 */
int main(int argc, char* argv[])
{
	for (int i = 0; i < argc; i++)
	{
		if (strcmp(argv[i], "-h") == 0 ||
			strcmp(argv[i], "-help") == 0)
		{
			PrintHelp();
			return 0;
		}

		if (strcmp(argv[i], "-i") == 0 ||
			strcmp(argv[i], "-intrinsic") == 0)
		{
			printIntrinsic = true;
		}

		if (strcmp(argv[i], "-d") == 0 ||
			strcmp(argv[i], "-debug") == 0)
		{
			printDebugPoints = true;
		}

		if (strcmp(argv[i], "-c") == 0 ||
			strcmp(argv[i], "-create") == 0)
		{
			CreateConfigFile();
			return 0;
		}
	}

	Settings s;

	const string inputSettingsFile = argc > 1 ? argv[1] : "default.xml";
	FileStorage fs(inputSettingsFile, FileStorage::READ); // Read the settings

	if (!fs.isOpened())
	{
		cout << "Could not open the configuration file: \"" << inputSettingsFile << "\"" << endl;
		return -1;
	}

	fs["Settings"] >> s;
	fs.release(); // close Settings file

	if (!s.goodInput)
	{
		cout << "Invalid input detected. Application stopping. " << endl;
		return -1;
	}

	vector<vector<Point2f> > imagePoints;
	Mat cameraMatrix, distCoeffs;
	Size imageSize;
	int mode = s.inputType == Settings::IMAGE_LIST ? CAPTURING : DETECTION;
	clock_t prevTimestamp = 0;
	const Scalar RED(0, 0, 255), GREEN(0, 255, 0);
	const char ESC_KEY = 27;
	vector<Point2f> pointBuf;
	bool blinkOutput, patternFound;

	while (true)
	{
		Mat view;
		pointBuf.clear();
		blinkOutput = false;
		patternFound = false;

		view = s.nextImage();

		if (view.empty()) // If no more images then run calibration, save and stop loop.
		{
			if (imagePoints.size() > 0)
				runCalibrationAndSave(s, imageSize, cameraMatrix, distCoeffs, imagePoints);
			break;
		}


		imageSize = view.size();  // Format input image.

		switch (s.calibrationPattern) // Find feature points on the input format
		{
		case Settings::CHESSBOARD:
			patternFound = findChessboardCorners(view, s.boardSize, pointBuf,
				CV_CALIB_CB_ADAPTIVE_THRESH | CV_CALIB_CB_FAST_CHECK | CV_CALIB_CB_NORMALIZE_IMAGE | CV_CALIB_CB_FILTER_QUADS);
			break;
		case Settings::CIRCLES_GRID:
			patternFound = findCirclesGrid(view, s.boardSize, pointBuf);
			break;
		case Settings::ASYMMETRIC_CIRCLES_GRID:
			patternFound = findCirclesGrid(view, s.boardSize, pointBuf, CALIB_CB_ASYMMETRIC_GRID);
			break;
		default:
			patternFound = false;
			break;
		}

		if (patternFound) // If done with success,
		{
			// improve the found corners' coordinate accuracy for chessboard
			if (s.calibrationPattern == Settings::CHESSBOARD)
			{
				Mat viewGray;
				cvtColor(view, viewGray, COLOR_BGR2GRAY);
				cornerSubPix(viewGray, pointBuf, Size(11, 11),
					Size(-1, -1), TermCriteria(CV_TERMCRIT_EPS + CV_TERMCRIT_ITER, 30, 0.1));
			}

			imagePoints.push_back(pointBuf);

			// Draw the corners.
			drawChessboardCorners(view, s.boardSize, Mat(pointBuf), patternFound);
			//imshow("Image View", view);
		}
	}

	// -----------------------Show the undistorted image for the image list ------------------------
	if (s.inputType == Settings::IMAGE_LIST)
	{
		Mat view, rview, map1, map2;
		vector<int> imgOutParams;
		imgOutParams.push_back(CV_IMWRITE_JPEG_QUALITY);
		imgOutParams.push_back(95);

		initUndistortRectifyMap(cameraMatrix, distCoeffs, Mat(),
			getOptimalNewCameraMatrix(cameraMatrix, distCoeffs, imageSize, 1, imageSize, 0),
			imageSize, CV_16SC2, map1, map2);

		String path = "";
		size_t loc = s.imageList[0].find_last_of("/\\");
		
		if (loc != string::npos)
			path = s.imageList[0].substr(0, s.imageList[0].find_last_of("/\\")) + "\\";

		// check if input image is available
		view = imread(path + "input.jpg", 1);
		
		// if no input then use the first calibration image
		if (view.empty())
			view = imread(s.imageList[0], 1);
		
		if (!view.empty())
		{
			remap(view, rview, map1, map2, INTER_LINEAR);
			//imshow("Image View", rview);
			imwrite(path + "out.jpg", rview, imgOutParams);

			if (printDebugPoints)
			{
				cout << "Debug Info" << endl;
				cout << "Input Points" << endl;
				cout << format(imagePoints[0], "csv") << endl;

				findChessboardCorners(rview, s.boardSize, pointBuf,
					CV_CALIB_CB_ADAPTIVE_THRESH | CV_CALIB_CB_FAST_CHECK | CV_CALIB_CB_NORMALIZE_IMAGE | CV_CALIB_CB_FILTER_QUADS);

				Mat rviewGray;
				cvtColor(rview, rviewGray, COLOR_BGR2GRAY);
				cornerSubPix(rviewGray, pointBuf, Size(11, 11),
					Size(-1, -1), TermCriteria(CV_TERMCRIT_EPS + CV_TERMCRIT_ITER, 30, 0.1));

				cout << "Undistorted Points" << endl;
				cout << format(pointBuf, "csv") << endl;
			}
		}
	}

	return 0;
}