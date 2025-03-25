module myViewCam
export myView

import VideoIO, ImageView;

function myView() 
        camera = VideoIO.opencamera();
        buf = VideoIO.read(camera);
        guidict = ImageView.imshow(buf);
        while !eof(camera)
            VideoIO.read!(camera, buf);
            ImageView.imshow(guidict["gui"]["canvas"], buf);
            sleep(0.00001);
        end
    end
	close(camera)
end
