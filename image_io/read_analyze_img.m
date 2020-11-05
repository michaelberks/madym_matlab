function [img] = read_analyze_img(filename, hdr_data, output_type)
% Read image file of Mayo Analyze 7.5 data set. 
% 
% Based on Matlab function analyze75read. Reads image data from the IMG file of
% an Analyze 7.5 format data set (a pair of FILENAME.HDR and FILENAME.IMG
% files).  For single-frame images, I is an M-by-N array where M is the
% number of rows and N is the number of columns. For multi-dimensional
% images, I can be an M-by-N-by-O or M-by-N-by-O-by-P array where M is
% the number of rows, N is the number of columns, O is the number of
% slices per volume and P is the number of volumes or time points. The
% data type of I is consistent with the image data type specified in 
% the metadata obtained from the header file. 
% 
% Inputs:
%     filename : str default None
%         filename to analyze image/header pair to read. Should be extension
%         free, or either of the filename.hdr or filename.img pair. Must be
%         set if hdr_data not set
%     hdr_data : SimpleNamespace, default None,
%         analyze format header data structure. If None will be read from header
%         file, in which case filename must not be None
%     output_type : np.dtype, default np.float64,
%         Numpy datatype output image array converted to. Use None to leave unchanged
%         from type specified in header data
%     scale : float default 1.0,
%         Value by which output array is scaled by, if != 1.0, ouput will be divided
%         (NOT multplied) by scale
%     flip_y : bool default True,
%         If true, flips the output image about array axis 0 (vertical flip)
%     flip_x : bool default False,
%         If true, flips the output image about array axis 1 (horizontal flip)
% 
% Outputs:
%     img : np.array
%         Numpy array of image data in data type specified by output data. Will
%         be 2, 3 or 4D (n_y, n_x, n_z = 1, n_v = 1)
% 
% Class support
% -------------
% I can be logical, uint8, int16, int32, single, or double. Complex and
% RGB data types are not supported.
% 
% Based on Matlab function ANALYZE75READ, copyright 2005-2011 The MathWorks, Inc.
     
%Get image and header filename
if exist('filename', 'var') && ~isempty(filename)
   [hdr_dir, base_filename,~] = fileparts(filename);
    hdr_filename = fullfile(hdr_dir, [base_filename '.hdr']); 
end

%Check if given metadata, if not, need to load header
if ~exist('hdr_data', 'var') || isempty(hdr_data)
    % Will crash if no filename, that's desired behaviour
    hdr_data = read_analyze_hdr(hdr_filename);
end
[hdr_dir, base_filename,~] = fileparts(hdr_data.Filename);
img_filename = fullfile(hdr_dir, [base_filename '.img']);

% Unpack data from the bytes buffer into a numpy array of correct
% shape and data format, based on info in image header 
if ~isfield(hdr_data, 'Dimensions') || length(hdr_data.Dimensions) < 4     
    error('Incorrect header metadata structure, missing dimensions');
end

% Get byte-order
if ~ismember(hdr_data.ByteOrder,{'ieee-le','ieee-be'})
    error('Byte order %s not recognised, should be ieee-le or ieee-be',...
        hdr_data.ByteOrder);
end

%Read binary contents of header file into buffer
img_fid = fopen(img_filename, 'r', hdr_data.ByteOrder);
    
%We unpack differently for sparse and full, if sparse
%not in header (it won't be unless image written by us) assume full
if isfield(hdr_data, 'Sparse') 
    sparse = hdr_data.Sparse;
else
    sparse = false;
end

% Obtain precision_string for reading in the data in the right format.
convert_to_binary = false;
switch hdr_data.ImgDataType
    case 'DT_BINARY'
        input_type = 'uint8';
        element_sz = 5;
        convert_to_binary = true;
    case 'DT_UNSIGNED_CHAR'
        input_type = 'uint8';
        element_sz = 5;
    case 'DT_SIGNED_SHORT'
        input_type = 'int16';
        element_sz = 6;
    case 'DT_SIGNED_INT'
        input_type = 'int32';
        element_sz = 8;
    case 'DT_FLOAT'
        input_type = 'float32';
        element_sz = 8;
    case 'DT_DOUBLE'
        input_type = 'float64';
        element_sz = 12;
    otherwise
        error('Data type %s not supported',hdr_data.ImgDataType); 
end

%Check if given output type, if not set to input type
if ~exist('output_type', 'var') || isempty(output_type)
    output_type = input_type;
end
precision = [input_type '=>' output_type];

if sparse
    %Create correct shape container and put data into
    %indexed elements
    if strcmpi(output_type, 'float32')
        zeros_type = 'single';
    elseif strcmpi(output_type, 'float64')
        zeros_type = 'double';
    else
        zeros_type = output_type;
    end
    img = zeros(hdr_data.Dimensions(:)', zeros_type);

    %Need to work out how many elements we have
    n_idx = hdr_data.ImgFileSize / element_sz;

    %Size of file should divide by element size into an into
    %throw error if not
    if rem(n_idx, 1)
        error('Image filesize (%d) not divisble by %d element size',...
            hdr_data.ImgFileSize, element_sz);
    end

    %Can have no non-zero elements, in which case return
    if n_idx > 0
        %Unpack the data values then indices, adding 1 to the indices
        vals = fread(img_fid, n_idx, precision);
        idx = fread(img_fid, n_idx, 'int32=>int32') + 1; 
        
        % As data stored in fortran form, this should just work
        img(idx) = vals;
    end

else
    % Compute number of elements to be read.
    count = prod(hdr_data.Dimensions);

    % Call struct unpack with the correct precision string
    img_data = fread(img_fid, count, precision);

    % As data stored in fortran form, this should just work
    img = reshape(img_data, hdr_data.Dimensions(:)');
end
fclose(img_fid);

%For binary data convert to bool
if convert_to_binary
    img = img ~= 0;
end

return