function save_vid(name,F)
	writerObj = VideoWriter(name);
	writerObj.FrameRate = 10;

	 % open the video writer
	open(writerObj);
	% write the frames to the video
	for i=1:length(F)
	    % convert the image to a frame
	    frame = F(i) ;    
	    writeVideo(writerObj, frame);
	end
	% close the writer object
	close(writerObj);
end