











#include "texture.h"

void ltexinit(texdata *lt)
{
	lt->data = NULL;
}

void ltexfree(texdata *lt)
{
	free(lt->data);
	lt->data = NULL;
}

static FILE* savepngfp = NULL;

void png_file_write(png_structp png_ptr, png_bytep data, png_size_t length) 
{
	fwrite(data, 1, length, savepngfp);
}

void png_file_flush(png_structp png_ptr)
{
	fflush(savepngfp);
}

int savepng(const char* full, unsigned char* data, int sizex, int sizey, int channels)
{
	//FILE *fp;
	png_structp png_ptr;
	png_infop info_ptr;
	//png_colorp palette;
	int color_type;
	int bit_depth;
	png_colorp palette;	
	png_color_16 Colors[PNG_MAX_PALETTE_LENGTH];
	png_byte Trans[PNG_MAX_PALETTE_LENGTH];
	int npal;
	png_bytep* row_pointers;
	int y;
	png_text text_ptr[1];
	char srcstr[123];

	/* Open the file */
	savepngfp = fopen(full, "wb");
	if (savepngfp == NULL)
		return 0;

	/* Create and initialize the png_struct with the desired error handler
	* functions.  If you want to use the default stderr and longjump method,
	* you can supply NULL for the last three parameters.  We also check that
	* the library version is compatible with the one used at compile time,
	* in case we are using dynamically linked libraries.  REQUIRED.
	*/
	png_ptr = png_create_write_struct(PNG_LIBPNG_VER_STRING,
		(png_voidp) NULL, NULL, NULL);

	if (png_ptr == NULL)
	{
		fclose(savepngfp);
		return 0;
	}

	/* Allocate/initialize the image information data.  REQUIRED */
	info_ptr = png_create_info_struct(png_ptr);
	if (info_ptr == NULL)
	{
		fclose(savepngfp);
		png_destroy_write_struct(&png_ptr,  NULL);
		return 0;
	}

	color_type = PNG_COLOR_TYPE_RGB;

	if(channels == 4)
		color_type = PNG_COLOR_TYPE_RGBA;
	if(channels == 1)
		color_type = PNG_COLOR_TYPE_GRAY;

	bit_depth = 8;

	png_set_IHDR(png_ptr, info_ptr, sizex, sizey, bit_depth, color_type, PNG_INTERLACE_NONE, PNG_COMPRESSION_TYPE_DEFAULT, PNG_FILTER_TYPE_DEFAULT);

	/* Set error handling.  REQUIRED if you aren't supplying your own
	* error handling functions in the png_create_write_struct() call.
	*/
	if (setjmp(png_jmpbuf(png_ptr)))
	{
		/* If we get here, we had a problem writing the file */
		fclose(savepngfp);
		png_destroy_write_struct(&png_ptr, &info_ptr);
		return 0;
	}

	/* One of the following I/O initialization functions is REQUIRED */

	/* Set up the output control if you are using standard C streams */
	//png_init_io(png_ptr, savepngfp);
	png_set_write_fn(png_ptr, NULL, png_file_write, png_file_flush);

	row_pointers = (png_bytep*) malloc(sizeof(png_bytep) * sizey);

	if(!row_pointers)
	{
//		OUTMEM();
		return 0;
	}

	for (y=0; y<sizey; y++)
		row_pointers[y] = (png_byte*)&data[y*sizex*channels];

	sprintf(srcstr, "Rendered using %s Version %d", "chart", -1);
#ifdef USESTEAM
	strcat(srcstr, " Authorized Steam Build");
#endif
	text_ptr[0].key = "Source";
	text_ptr[0].text = srcstr;
	text_ptr[0].compression = PNG_TEXT_COMPRESSION_NONE;
	// text_ptr[2].compression = PNG_TEXT_COMPRESSION_zTXt;
#ifdef PNG_iTXt_SUPPORTED
	text_ptr[0].lang = NULL;
#endif
	text_ptr[0].lang_key = NULL;

	png_set_text(png_ptr, info_ptr, text_ptr, 1);

	png_write_info(png_ptr, info_ptr);

	png_write_image(png_ptr, row_pointers);
	png_write_end(png_ptr, NULL);

	//for (y=0; y<sizey; y++)
	//   free(row_pointers[y]);
	free(row_pointers);


	/* This is the easy way.  Use it if you already have all the
	* image info living in the structure.  You could "|" many
	* PNG_TRANSFORM flags into the png_transforms integer here.
	*/
	//png_write_png(png_ptr, info_ptr, NULL, NULL);

	/* If you png_malloced a palette, free it here (don't free info_ptr->palette,
	* as recommended in versions 1.0.5m and earlier of this example; if
	* libpng mallocs info_ptr->palette, libpng will free it).  If you
	* allocated it with malloc() instead of png_malloc(), use free() instead
	* of png_free().
	*/
	//png_free(png_ptr, palette);
	//palette = NULL;

	/* Similarly, if you png_malloced any data that you passed in with
	* png_set_something(), such as a hist or trans array, free it here,
	* when you can be sure that libpng is through with it.
	*/
	//png_free(png_ptr, trans);
	//trans = NULL;
	/* Whenever you use png_free() it is a good idea to set the pointer to
	* NULL in case your application inadvertently tries to png_free() it
	* again.  When png_free() sees a NULL it returns without action, thus
	* avoiding the double-free security problem.
	*/

	/* Clean up after the write, and free any memory allocated */
	png_destroy_write_struct(&png_ptr, &info_ptr);

	/* Close the file */
	fclose(savepngfp);

	/* That's it */
	return 1;
}

static FILE* loadpngfp = NULL;

void png_file_read(png_structp png_ptr, png_bytep data, png_size_t length)
{
	fread(data, 1, length, loadpngfp);
}

texdata *loadpng(const char *file)
{
	texdata *ltex = NULL;
	//header for testing if it is a png
	png_byte header[8];
	int is_png;
	png_structp png_ptr;
	png_infop info_ptr;
	png_infop end_info;
	int i;

	//variables to pass to get info
	int bit_depth, color_type;
	png_uint_32 twidth, theight;
	png_bytep trans_alpha = NULL;
	int num_trans = 0;
	png_color_16p trans_color = NULL;
	int row_bytes;
	png_bytep *row_pointers;

	loadpngfp = fopen(file, "rb");

	if (!loadpngfp)
	{
		return NULL;
	}

	fread(header, 8, 1, loadpngfp);
	//g_src.read((void*)header, 8);

	//test if png
	is_png = !png_sig_cmp(header, 0, 8);
	if (!is_png)
	{
	//	fprintf(g_applog,
	//		"Not a png file : %s %d,%d,%d,%d,%d,%d,%d,%d \r\n", file,
	//		(int)header[0], (int)header[1], (int)header[2], (int)header[3],
	//		(int)header[4], (int)header[5], (int)header[6], (int)header[7]);
		fclose(loadpngfp);
		return NULL;
	}

	//create png struct
	png_ptr = png_create_read_struct(PNG_LIBPNG_VER_STRING, NULL, NULL, NULL);
	if (!png_ptr)
	{
		//fprintf(g_applog, "Unable to create png struct : %s \r\n", file);
		fclose(loadpngfp);
		return NULL;
	}

	//create png info struct
	info_ptr = png_create_info_struct(png_ptr);
	if (!info_ptr)
	{
		png_destroy_read_struct(&png_ptr, (png_infopp)NULL, (png_infopp)NULL);
		//fprintf(g_applog, "Unable to create png info : %s \r\n", file);
		fclose(loadpngfp);
		return NULL;
	}

	//create png info struct
	end_info = png_create_info_struct(png_ptr);
	if (!end_info)
	{
		png_destroy_read_struct(&png_ptr, &info_ptr, (png_infopp)NULL);
	//	fprintf(g_applog, "Unable to create png end info : %s \r\n", file);
		fclose(loadpngfp);
		return NULL;
	}

	//png error stuff, not sure libpng man suggests this.
	if (setjmp(png_jmpbuf(png_ptr)))
	{
	//	fprintf(g_applog, "Error during setjmp : %s \r\n", file);
		png_destroy_read_struct(&png_ptr, &info_ptr, &end_info);
		fclose(loadpngfp);
		return NULL;
	}

	png_set_read_fn(png_ptr, NULL, png_file_read);

	//let libpng know you already read the first 8 bytes
	png_set_sig_bytes(png_ptr, 8);

	// read all the info up to the image data
	png_read_info(png_ptr, info_ptr);

	// get info about png
	png_get_IHDR(png_ptr, info_ptr, &twidth, &theight, &bit_depth, &color_type, NULL, NULL, NULL);

	//	dbool alphaFlag;

	ltex = (texdata*)malloc(sizeof(texdata));
	ltexinit(ltex);

	ltex->sizex = twidth;
	ltex->sizey = theight;

	switch (color_type)
	{
	case PNG_COLOR_TYPE_RGBA:
		ltex->channels = 4;
		break;
	case PNG_COLOR_TYPE_RGB:
		ltex->channels = 3;
		break;
	case PNG_COLOR_TYPE_PALETTE:
	{
		png_set_palette_to_rgb(png_ptr);
		ltex->channels = 3;
#if 0
		png_get_tRNS(png_ptr, info_ptr, &trans_alpha, &num_trans, &trans_color);
		if (trans_alpha != NULL)
			alphaFlag = ectrue;
		else
			alphaFlag = ecfalse;
		if (alphaFlag)
			ltex->channels = 4;
#endif

		//if(alphaFlag)
		if (png_get_valid(png_ptr, info_ptr, PNG_INFO_tRNS))
		{
			ltex->channels = 4;
			png_set_tRNS_to_alpha(png_ptr);
		}
	}
	break;
	default:
	//	fprintf(g_applog,
//			"%s color type %d not supported \r\n", file, (int)png_get_color_type(png_ptr, info_ptr));

		png_destroy_read_struct(&png_ptr, &info_ptr, NULL);
		ltexfree(ltex);
		free(ltex);
		fclose(loadpngfp);
		return NULL;
	}

	// Update the png info struct.
	png_read_update_info(png_ptr, info_ptr);

	// Row size in bytes.
	row_bytes = png_get_rowbytes(png_ptr, info_ptr);

	ltex->data = (unsigned char*)malloc(row_bytes * ltex->sizey);

	// Allocate the image_data as a big block, to be given to opengl
	if (!ltex->data)
	{
		//clean up memory and close stuff
		png_destroy_read_struct(&png_ptr, &info_ptr, &end_info);
//		fprintf(g_applog, "Unable to allocate image_data while loading %s \r\n", file);
		ltexfree(ltex);
		free(ltex);
		fclose(loadpngfp);
		return NULL;
	}

	//row_pointers is for pointing to image_data for reading the png with libpng
	row_pointers = (png_bytep*)malloc(sizeof(png_bytep)*ltex->sizey);
	if (!row_pointers)
	{
		//clean up memory and close stuff
		png_destroy_read_struct(&png_ptr, &info_ptr, &end_info);
		ltexfree(ltex);
		free(ltex->data);
	//	fprintf(g_applog, "Unable to allocate row_pointer while loading %s \r\n", file);
		fclose(loadpngfp);
		return NULL;
	}

	// set the individual row_pointers to point at the correct offsets of image_data
	for (i = 0; i < ltex->sizey; ++i)
		row_pointers[i] = ltex->data + i * row_bytes;

	//read the png into image_data through row_pointers
	png_read_image(png_ptr, row_pointers);

	//clean up memory and close stuff
	png_destroy_read_struct(&png_ptr, &info_ptr, &end_info);
	free(row_pointers);
	fclose(loadpngfp);

	return ltex;
}