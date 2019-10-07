
var glcanvas;
var gl;
var canvas;
var ctx;

var voxdata;

function init(){	
	glcanvas=document.getElementById("glcanvas");
	gl=glcanvas.getContext("webgl");
	
	mycanvas=document.getElementById("mycanvas");
	ctx=mycanvas.getContext("2d");
	
	
	//voxel data. maybe better to have octree or similar. for now, use array of arrays. implementation maybe sparse. maybe good to use bits in a number (32 bit?)
	var blocksize = 64;	//64x64x64 = 512x512, which seems manageable to draw slices onto canvas
	var canvassize = 512;
	
	voxdata = [];
	//generate arrays 
	for (var ii=0;ii<blocksize;ii++){
		var slicedata=[];
		voxdata.push(slicedata);
		for (var jj=0;jj<blocksize;jj++){
			slicedata.push([]);	//could push new Array(blocksize). use empty array if fill with sparse data. unimportant for small blocksize 
		}
	}
	
	//fill with data
	function makeVoxdataForFunc(thefunction){
		for (var ii=0;ii<blocksize;ii++){
			var slicedata = voxdata[ii];
			for (var jj=0;jj<blocksize;jj++){
				var stripdata = slicedata[jj];
				for (var kk=0;kk<blocksize;kk++){
					stripdata[kk] = thefunction(ii,jj,kk);	//if know starting from empty array, can use push. if want sparse array, should use conditional. 
				}											//note some functions may benefit from currying eg functionForXY= thefunction(ii,jj), 
			}
		}
	}
	
	makeVoxdataForFunc(sinesfunction);
	
	function sinesfunction(ii,jj,kk){ 
		return Math.sin(ii)+Math.sin(jj)+Math.sin(kk) > 0;
	}
	
	console.log(voxdata);
	
	//draw to canvas to test algo.
	//start with a single slice
	var slicedata = voxdata[0];
	var imgData = ctx.createImageData(blocksize, blocksize);
	var imgDataData = imgData.data;
	var idx=0;
	for (var jj=0;jj<blocksize;jj++){
		var stripdata = slicedata[jj];
		for (var kk=0;kk<blocksize;kk++){
			var color = stripdata[kk] ? 255 : 0;
			imgDataData[idx]=color;
			imgDataData[idx+1]=color;
			imgDataData[idx+2]=color;
			imgDataData[idx+3]=255;
			idx+=4;
		}
	}
	ctx.putImageData(imgData, 0, 0);
}
