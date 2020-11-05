function hdr_data = read_analyze_hdr(filename, byte_order)
%read_analyze_hdr Read hdr_data from header file of Mayo Analyze 7.5 data set.
%
%   HDR_DATA = READ_ANALYZE_HDR(FILENAME) reads the HDR file of an Analyze 7.5
%   format data set (a pair of FILENAME.HDR and FILENAME.IMG files) and
%   returns a structure HDR_DATA whose fields contain information about the
%   data set. FILENAME is a string that specifies the name of the Analyze
%   file pair.
%
%   HDR_DATA = READ_ANALYZE_HDR(FILENAME, BYTE_ORDER) attempts to
%   read the HDR file with the byte-ordering specified in the string
%   ENDIAN. If the specified ENDIAN value results in a read error,
%   read_analyze_hdr generates a warning and attempts to read the HDR file
%   with the opposite ByteOrder format. Valid values for ENDIAN can be
%   'ieee-le' (little endian) or 'ieee-be' (big endian).
%
%   Example 1
%   ---------
%   Use read_analyze_hdr to read the header file brainMRI.hdr. 
%
%       hdr_data = read_analyze_hdr('brainMRI.hdr');
%
%   Example 2
%   ---------
%   Specify the Byte ordering of the data set. Attempt reading the
%   hdr_data from the header file using the specified Byte ordering.
%
%       hdr_data = read_analyze_hdr('brainMRI', 'ieee-le');
%
%   See also READ_ANALYZE_IMG.

%   Following is a partial list of fields in the HDR_DATA structure:
%
%       Filename          A string containing the name of the file
%
%       FileModDate       A string containing the modification date of
%                         the file
%
%       HdrFileSize       An integer indicating the size of the HDR file in
%                         bytes
%
%       ImgFileSize       An integer indicating the size of the IMG file in
%                         bytes
%
%       Format            A string containing the file format. This value is
%                         set to 'Analyze' for valid Analyze data sets
%
%       FormatVersion     A string or number specifying the file format
%                         version
%
%       Width             An integer indicating the width of the image
%                         in pixels
%
%       Height            An integer indicating the height of the image
%                         in pixels
%
%       BitDepth          An integer indicating the number of bits per
%                         pixel
%
%       ColorType         A string indicating the type of image; either
%                         'truecolor' for a truecolor (RGB) image, or
%                         'grayscale' for a grayscale image,
%
%       ByteOrder         A string containing the byte-ordering used to
%                         successfully read in the HDR file.
%
%       HdrDataType       Data type of the HDR file.
%
%       DatabaseName      Name of the image database.
%
%       Extents           An integer which is a required field in the header
%                         file. This value should be 16384.
%
%       SessionError      An integer indicating session error number.
%
%       Regular           A character indicating whether or not all images
%                         and volumes are of the same size. A value '1'
%                         indicates that the data is regular while '0'
%                         indicates the data is not regular.
%
%       Dimensions        A vector providing information on the image
%                         dimensions. The vector is of the form
%                             [X Y Z T]
%                         X gives the X dimension of the image, i.e. the
%                         number of pixels in an image row.
%                         Y gives the Y dimension of the image, i.e. the
%                         number of pixels in an image column.
%                         Z gives the volume Z dimension, i.e. the number of
%                         slices in a volume.
%                         T indicates the time points, i.e. the number of
%                         volumes in the dataset.
%                         Dimensions vector only returns non-zero entries.
%
%       VoxelUnits        Spatial units of measure for a voxel.
%
%       CalibrationUnits  Name of the calibration unit.
%
%       ImgDataType       Data type of the IMG file.
%
%       PixelDimensions   A vector providing information on the pixel
%                         dimensions. PixelDimensions is parallel to the
%                         Dimensions field, providing real world
%                         measurements in mm. The vector is of the form
%                             [Xp Yp Zp Tp]
%                         Xp provides the voxel width in mm.
%                         Yp provides the voxel height in mm.
%                         Zp provides the slice thickness in mm.
%                         Tp provides the time points in ms.
%                         PixelDimensions vector only returns non-zero
%                         entries.
%
%       VoxelOffset       The byte offset in the image file at which voxels
%                         start. This value may be negative to specify that
%                         the absolute value is applied for every image in
%                         the file.
%
%       CalibrationMax    Maximum Calibration value.
%
%       CalibrationMin    Minimum Calibration value.
%
%       GlobalMax         Global Maximum. The maximum pixel values for the
%                         entire dataset.
%
%       GlobalMin         Global Minimum. The minimum pixel values for the
%                         entire dataset.
%
%       Descriptor        Data description.
%
%       Orientation       Slice orientation for the dataset.
%

%   Copyright 2005-2011 The MathWorks, Inc.

%Add hdr extension (in case .img file given as input)
[hdr_dir, base_filename,~] = fileparts(filename);
hdr_filename = fullfile(hdr_dir, [base_filename '.hdr']);
if ~exist(hdr_filename, 'file')
    error('%s does not exist', hdr_filename);
end

% Check byte order value
user_supplied = false;
if exist('byte_order', 'var') && ~isempty(byte_order)
    user_supplied = true;

    if ~ismember(byte_order, {'ieee-le', 'ieee-be'})
        error([byte_order ' not recognised, must be ieee-le or ieee-be']);
    end
end  

% Remember if we have already warned the user about truncated header file.
already_warned = false;

%--------------------------------------------------------------------------------
%Set up hdr_data
hdr_data.Filename = hdr_filename;
d = dir(hdr_filename);
hdr_data.FileModDate = d.date;

% Image File Size will be obtained below after constructing the
% Image filename from the header filename.
hdr_data.Format = 'Analyze';
hdr_data.FormatVersion = '7.5';
hdr_data.Width = [];
hdr_data.Height = [];
hdr_data.BitDepth = [];
hdr_data.ColorType = 'unknown';

% Construct Image filename using Filename obtained from
% Metadata struct.
img_filename = fullfile(hdr_dir, [base_filename '.img']);  

% Obtain Image file size.
try
    d = dir(img_filename);
    hdr_data.ImgFileSize = d.bytes;
catch
    hdr_data.ImgFileSize = [];
    warning('%s missing image part %s', hdr_filename, img_filename);
end

%-------------------------------------------------------------------------
%Check whether buffer is little or big-endian
%Read headerSize with both little and big endian
%Read binary contents of header file into buffer
fid_le = fopen(hdr_filename, 'r', 'l');
fid_be = fopen(hdr_filename, 'r', 'b');

headerSize_le = fread(fid_le, 1, 'int32=>int32');
headerSize_be = fread(fid_be, 1, 'int32=>int32');

% Possible exted header size
min_size = 348; 
max_size = 2000;

% headerSize should be within the extedRange. Use that to check
% if incorrect ByteOrder was used to open and read the file.
if (min_size <= headerSize_le) && (headerSize_le <= max_size)
    fid = fid_le;
    fclose(fid_be);
    hdr_data.HdrFileSize = headerSize_le;
    hdr_data.ByteOrder = 'ieee-le';

elseif (min_size <= headerSize_be) && (headerSize_be <= max_size)
    fid = fid_be;
    fclose(fid_le);
    hdr_data.HdrFileSize = headerSize_be;
    hdr_data.ByteOrder = 'ieee-be';

else        
    % We have tried reading the file with both ByteOrder
    % formats. Generate error
    fclose(fid_le);
    fclose(fid_be);
    error('Analyze header size wrong, is this really an hdr file')
end

% Generate warning if incorrect ByteOrder was provided by user.
if user_supplied && ~strcmpi(hdr_data.ByteOrder, byte_order)
    warning('User supplied %s was incorrect, using %s.',...
        hdr_data.ByteOrder, byte_order);
end
%------------------------------------------------------------------------------------
% Read the HeaderKey information from HDR file.
% Read all information in the HeaderKey structure.
hdr_data.HdrDataType  = unpackVerified(10, 's');
hdr_data.DatabaseName = unpackVerified(18, 's');
hdr_data.Extents      = unpackVerified(1,  'i');
hdr_data.SessionError = unpackVerified(1,  'h');
hdr_data.Regular      = unpackVerified(1,  's')  == 'r';

%Advance one position for an unused character.
fseek(fid,1,'cof');

% Read the ImgDimension information from HDR file.
dims = unpackVerified(8, 'h');

% Return useful dimension information.
dims(1) = [];
hdr_data.Dimensions = dims(dims ~= 0);
hdr_data.Width = hdr_data.Dimensions(1);
hdr_data.Height= hdr_data.Dimensions(2);
hdr_data.VoxelUnits  = unpackVerified(4, 's');
hdr_data.CalibrationUnits = unpackVerified(8, 's');

% Advance 2 positions for an unused field.
fseek(fid,2,'cof');

% Parse image data type
ImgDataType = unpackVerified(1, 'h');

%Check if sparse
%Our sparse format adds 5 to the supported data types,
%(this means we can sparsify masks, if we added 1, then
% a sparse BINARY type would have code 2 == UNSIGNED_CHAR)
hdr_data.Sparse = false;
if ImgDataType == 6 || rem(ImgDataType, 2)
    hdr_data.Sparse = true;
    ImgDataType = ImgDataType - 5;
end

switch ImgDataType
    case int16(0)
        hdr_data.ImgDataType = 'DT_UNKNOWN';
    case int16(1)
        hdr_data.ImgDataType = 'DT_BINARY';
        hdr_data.ColorType = 'grayscale';
    case int16(2)
        hdr_data.ImgDataType = 'DT_UNSIGNED_CHAR';
        hdr_data.ColorType = 'grayscale';
    case int16(4)
        hdr_data.ImgDataType = 'DT_SIGNED_SHORT';
        hdr_data.ColorType = 'grayscale';
    case int16(8)
        hdr_data.ImgDataType = 'DT_SIGNED_INT';
        hdr_data.ColorType = 'grayscale';
    case int16(16)
        hdr_data.ImgDataType = 'DT_FLOAT';
        hdr_data.ColorType = 'grayscale';
    case int16(32)
        hdr_data.ImgDataType = 'DT_COMPLEX';
        hdr_data.ColorType = 'grayscale';
    case int16(64)
        hdr_data.ImgDataType = 'DT_DOUBLE';
        hdr_data.ColorType = 'grayscale';
    case int16(128)
        hdr_data.ImgDataType = 'DT_RGB';
        hdr_data.ColorType = 'truecolor';
    case int16(255)
        hdr_data.ImgDataType = 'DT_ALL';
end  % switch

%Read bit depth
hdr_data.BitDepth   = unpackVerified(1, 'h');

%Advance 2 positions for an unused field.
fseek(fid,2,'cof');
pix_dims   = unpackVerified(8, 'f');
hdr_data.PixelDimensions = pix_dims(pix_dims ~= 0);
hdr_data.VoxelOffset     = unpackVerified(1, 'f');

%Advance 12 positions for an unused field.
fseek(fid,12,'cof');
hdr_data.CalibrationMax = unpackVerified(1, 'f');
hdr_data.CalibrationMin = unpackVerified(1, 'f');
hdr_data.Compressed     = unpackVerified(1, 'f');
hdr_data.Verified       = unpackVerified(1, 'f');
hdr_data.GlobalMax      = unpackVerified(1, 'i');
hdr_data.GlobalMin      = unpackVerified(1, 'i') ;

% Read the DataHistory information from HDR file.
hdr_data.Descriptor   = unpackVerified(80, 's');
hdr_data.AuxFile      = unpackVerified(24, 's');
Orientation           = unpackVerified(1, 'B');

switch Orientation
    case 0
        hdr_data.Orientation = 'Transverse unflipped';
    case 1
        hdr_data.Orientation = 'Coronal unflipped';
    case 2
        hdr_data.Orientation = 'Sagittal unflipped';
    case 3
        hdr_data.Orientation = 'Transverse flipped';
    case 4
        hdr_data.Orientation = 'Coronal flipped';
    case 5
        hdr_data.Orientation = 'Sagittal flipped';
    otherwise
        hdr_data.Orientation = 'Orientation unavailable';
end

% Various scanner generated fields
hdr_data.Originator   = unpackVerified(10, 's');
hdr_data.Generated    = unpackVerified(10, 's');
hdr_data.Scannumber   = unpackVerified(10, 's');
hdr_data.PatientID    = unpackVerified(10, 's');
hdr_data.ExposureDate = unpackVerified(10, 's');
hdr_data.ExposureTime = unpackVerified(10, 's');

%Advance 3 positions for an unused field.
fseek(fid,3,'cof');
hdr_data.Views          = unpackVerified(1, 'i');
hdr_data.VolumesAdded   = unpackVerified(1, 'i');
hdr_data.StartField     = unpackVerified(1, 'i');
hdr_data.FieldSkip      = unpackVerified(1, 'i');
hdr_data.OMax           = unpackVerified(1, 'i');
hdr_data.OMin           = unpackVerified(1, 'i');
hdr_data.SMax           = unpackVerified(1, 'i');
hdr_data.SMin           = unpackVerified(1, 'i');
fclose(fid);

%We're done, return header data
return

%%
%------------------------------------------------------------------------
%%%
%%% Define nested helper function unpackVerified
%%%
    function [out] = unpackVerified(count, precision)
        % This function reads the specified number of bytes using unpack and
        % checks for premature EOF. In that case, a warning is generated the
        % first time this is encountered in the file.
        switch precision
            case 's' %string
                precision_str = 'uchar=>char';
            case 'B' %uint8
                precision_str = 'uint8=>uint8';
            case 'h' %short
                precision_str = 'int16=>int16';
            case 'i' %int
                precision_str = 'int32=>int32';
            case 'f' %float
                precision_str = 'float32=>float32';
            otherwise
                error('Precision %s not recognised, should be s, B, h, i or f',...
                    precision);
        end
        temp = fread(fid, count, precision_str);
        if isempty(temp)
            if ~already_warned
                warning('Truncated header file');

                % Set alreadyWarned to true so that we don't warn again.
                already_warned = true;
            end  % if

        end  % if
        
        if precision == 's'
            out = deblank(temp);
        else
            out = temp;
        end
    end
%%%
end